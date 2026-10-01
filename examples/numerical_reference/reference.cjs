'use strict';
// Generic, nondimensional numerical reference. Not a host or pathogen model.
const crypto = require('node:crypto');
const ENGINE = {id:'pocket-numerical-reference',version:'1.0.0',scope:'synthetic fields only; no biological calibration'};
const digest = x => crypto.createHash('sha256').update(x).digest('hex');
function randomFor(seed,id,patch){
  let v=parseInt(digest(JSON.stringify([seed,id,patch])).slice(0,8),16);
  return ()=>{v=(Math.imul(v,1664525)+1013904223)>>>0;return v/4294967296;};
}
function finite(v,label,min=0){if(!Number.isFinite(v)||v<min)throw Error('Invalid '+label);return v;}
function configuration(input={}){
  const c={nx:12,ny:24,width:1,depth:1,dt:.01,duration:.4,seed:19,diffusion:.002,resourceDiffusion:.01,
    resourceInitial:[.1,.08],boundaryRate:0,boundaryReservoir:[.2,.1],
    taxa:[{id:'synthetic-a',rate:.7,mortality:.04,initial:.01,costs:[.5,0]},
          {id:'synthetic-b',rate:.6,mortality:.04,initial:.012,costs:[.2,.3],supportFrom:'synthetic-a',supportScale:.2}],...input};
  for(const k of ['nx','ny','seed'])if(!Number.isSafeInteger(c[k]))throw Error('Integer '+k+' required');
  if(c.nx<2||c.ny<2||c.nx*c.ny>100000)throw Error('Invalid grid');
  for(const k of ['width','depth','dt'])finite(c[k],k,Number.MIN_VALUE);
  for(const k of ['duration','diffusion','resourceDiffusion','boundaryRate'])finite(c[k],k);
  if(c.duration/c.dt>200000||c.duration/c.dt>Number.MAX_SAFE_INTEGER)throw Error('Run too large');
  for(const key of ['resourceInitial','boundaryReservoir']){
    if(!Array.isArray(c[key])||c[key].length!==2)throw Error('Two resources required');c[key].forEach(v=>finite(v,key));
  }
  if(!Array.isArray(c.taxa)||!c.taxa.length||c.taxa.length>10)throw Error('Invalid synthetic panel');
  c.taxa=c.taxa.map(t=>({...t,costs:[...t.costs]})).sort((a,b)=>a.id.localeCompare(b.id));
  const ids=new Set();
  for(const t of c.taxa){
    if(!/^synthetic-[a-z0-9-]+$/.test(t.id)||ids.has(t.id))throw Error('Unique synthetic IDs required');ids.add(t.id);
    for(const k of ['rate','mortality','initial'])finite(t[k],k);
    if(t.costs.length!==2||!t.costs.some(v=>v>0))throw Error('Resource cost required');t.costs.forEach(v=>finite(v,'cost'));
    finite(t.supportScale??0,'supportScale');
  }
  for(const t of c.taxa)if(t.supportFrom&&!ids.has(t.supportFrom))throw Error('Unknown synthetic support ID');
  const maxDiff=Math.max(c.diffusion,c.resourceDiffusion),outflow=2*maxDiff*((c.nx/c.width)**2+(c.ny/c.depth)**2);
  if(c.dt*outflow>.9)throw Error('Combined diffusion positivity bound exceeded; reduce dt');
  return c;
}
function integral(a,c){return a.reduce((s,v)=>s+v,0)*c.width*c.depth/(c.nx*c.ny);}
function diffuse(a,c,D,dt){
  const out=Float64Array.from(a),fx=D*dt*(c.nx/c.width)**2,fy=D*dt*(c.ny/c.depth)**2;
  for(let y=0;y<c.ny;y++)for(let x=0;x<c.nx;x++){
    const i=y*c.nx+x;
    if(x+1<c.nx){const q=fx*(a[i+1]-a[i]);out[i]+=q;out[i+1]-=q;}
    if(y+1<c.ny){const q=fy*(a[i+c.nx]-a[i]);out[i]+=q;out[i+c.nx]-=q;}
  }
  return out;
}
function initialize(c){
  const n=c.nx*c.ny;
  const biomass=Object.fromEntries(c.taxa.map(t=>{
    const patches=Array.from({length:7},(_,p)=>{const r=randomFor(c.seed,t.id,'patch-'+p);return [r(),r(),.06+.02*r()];});
    const a=Float64Array.from({length:n},(_,i)=>{
      const x=(i%c.nx+.5)/c.nx,y=(Math.floor(i/c.nx)+.5)/c.ny;
      if(c.smoothInitial)return t.initial*(1+.2*Math.cos(Math.PI*x)*Math.cos(Math.PI*y));
      return t.initial*(.25+patches.reduce((s,[u,v,r])=>s+Math.exp(-((x-u)**2+(y-v)**2)/(2*r*r)),0));
    });return [t.id,a];
  }));
  return {biomass,resources:c.resourceInitial.map(v=>new Float64Array(n).fill(v))};
}
function step(state,c,dt=c.dt){
  finite(dt,'step',Number.MIN_VALUE);if(dt>c.dt)throw Error('Step exceeds validated dt');
  const b=Object.fromEntries(c.taxa.map(t=>[t.id,diffuse(state.biomass[t.id],c,c.diffusion,dt)]));
  const r=state.resources.map(a=>diffuse(a,c,c.resourceDiffusion,dt));
  const cellArea=c.width*c.depth/(c.nx*c.ny);
  const ledger={birth:0,loss:0,consumed:[0,0],exchange:[0,0],maximumResourceOverdraw:0};
  // Exact nonnegative relaxation on top row; exchange is measured from actual delta.
  const fraction=-Math.expm1(-c.boundaryRate*dt);
  for(let k=0;k<2;k++)for(let i=0;i<c.nx;i++){
    const delta=fraction*(c.boundaryReservoir[k]-r[k][i]);r[k][i]+=delta;ledger.exchange[k]+=delta*cellArea;
  }
  for(let i=0;i<r[0].length;i++){
    // All rates read this same state. No newly credited birth supports another field.
    const old=c.taxa.map(t=>b[t.id][i]);
    const proposed=c.taxa.map((t,s)=>{
      const fit=Math.min(...t.costs.map((cost,k)=>cost>0?r[k][i]/(.1+r[k][i]):1));
      const donor=t.supportFrom?old[c.taxa.findIndex(v=>v.id===t.supportFrom)]:0;
      return dt*t.rate*old[s]*fit*(1+(t.supportScale??0)*donor/(.1+donor));
    });
    const demand=[0,0];for(let s=0;s<c.taxa.length;s++)for(let k=0;k<2;k++)demand[k]+=proposed[s]*c.taxa[s].costs[k];
    // Common resource-allocation factor: conservative and symmetric, not maximally efficient.
    // A 16-epsilon reserve prevents subtraction roundoff; it is never credited as growth.
    let scale=1;for(let k=0;k<2;k++)if(demand[k]>0)scale=Math.min(scale,r[k][i]*(1-16*Number.EPSILON)/demand[k]);
    const consumed=[0,0];
    for(let s=0;s<c.taxa.length;s++){
      const t=c.taxa[s],birth=proposed[s]*scale,loss=old[s]*(-Math.expm1(-t.mortality*dt));
      b[t.id][i]=old[s]+birth-loss;ledger.birth+=birth*cellArea;ledger.loss+=loss*cellArea;
      for(let k=0;k<2;k++)consumed[k]+=birth*t.costs[k];
    }
    for(let k=0;k<2;k++){ledger.maximumResourceOverdraw=Math.max(ledger.maximumResourceOverdraw,consumed[k]-r[k][i]);r[k][i]-=consumed[k];ledger.consumed[k]+=consumed[k]*cellArea;}
  }
  for(const a of [...Object.values(b),...r])for(const v of a)if(!Number.isFinite(v)||v<0)throw Error('Nonfinite/negative state; no silent clipping');
  return {state:{biomass:b,resources:r},ledger};
}
function simulate(options={}){
  const c=configuration(options);let state=initialize(c);
  const mass=s=>({biomass:Object.values(s.biomass).reduce((v,a)=>v+integral(a,c),0),resources:s.resources.map(a=>integral(a,c))});
  const initial=mass(state),ledger={birth:0,loss:0,consumed:[0,0],exchange:[0,0],maximumResourceOverdraw:0};
  const steps=Math.ceil(c.duration/c.dt-1e-12);
  for(let j=0;j<steps;j++){
    const result=step(state,c,Math.min(c.dt,c.duration-j*c.dt));state=result.state;
    ledger.birth+=result.ledger.birth;ledger.loss+=result.ledger.loss;
    ledger.maximumResourceOverdraw=Math.max(ledger.maximumResourceOverdraw,result.ledger.maximumResourceOverdraw);
    for(let k=0;k<2;k++){ledger.consumed[k]+=result.ledger.consumed[k];ledger.exchange[k]+=result.ledger.exchange[k];}
  }
  const final=mass(state),all=[...Object.values(state.biomass),...state.resources];
  const minimum=Math.min(...all.map(a=>Math.min(...a))),finiteState=all.every(a=>a.every(Number.isFinite));
  return {engine:ENGINE,config:c,state,initial,final,ledger,audit:{minimum,finite:finiteState,
    biomassResidual:final.biomass-initial.biomass-ledger.birth+ledger.loss,
    resourceResiduals:final.resources.map((v,k)=>v-initial.resources[k]-ledger.exchange[k]+ledger.consumed[k])}};
}
module.exports={ENGINE,digest,randomFor,configuration,integral,diffuse,initialize,step,simulate};
