import {test} from 'node:test';
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {readFile} from 'node:fs/promises';
import {buildCrevice,microbialObject,microbialPlacements,creviceGap,creviceLandmarks,creviceCamera,creviceProfiles,creviceLabels,crevicePoint,toothWall,CREVICE_CLAIM,CREVICE_SWEEP} from '../dist/pocket-core/pocket-crevice.mjs';
import {TAXA,PRESETS} from '../dist/pocket-core/pocket-science.mjs';
const hash=g=>createHash('sha256').update(Buffer.from(g.getAttribute('position').array.buffer)).digest('hex');
test('all nine objects are deterministic, distinct, finite and source identified',()=>{
 const hashes=[];
 for(const t of TAXA){const a=microbialObject(t.id),b=microbialObject(t.id);assert.equal(hash(a),hash(b));hashes.push(hash(a));for(const key of ['position','normal'])assert.ok(a.getAttribute(key).array.every(Number.isFinite));a.dispose();b.dispose();}
 assert.equal(new Set(hashes).size,9);assert.throws(()=>microbialObject('unresolved'),/Unknown microbial/);
 for(const health of ['healthy','periodontitis']){const parts=buildCrevice(health),cells=parts.filter(p=>p.taxon);assert.deepEqual(cells.map(c=>c.taxon),TAXA.map(t=>t.id));for(const c of cells)assert.match(c.source,/^https:\/\//);for(const p of parts){assert.ok(p.geometry.getAttribute('normal').array.every(Number.isFinite));p.geometry.dispose();}}
});
test('entrance stays narrower than the authored interior and landmarks remain separate',()=>{
 const original=JSON.stringify(PRESETS);
 for(const health of ['healthy','periodontitis']){const {margin,attachment,crest}=creviceLandmarks(health);assert.ok(margin>attachment&&attachment>crest);assert.ok(creviceGap(margin,health)<creviceGap((margin+attachment)/2,health));assert.ok(creviceGap(attachment,health)<creviceGap((margin+attachment)/2,health));for(const p of microbialPlacements(health))assert.ok(p.position[1]<margin&&p.position[1]>attachment);}
 assert.equal(JSON.stringify(PRESETS),original);
});
test('crevice cameras are distinct and the isolated form focus follows selected identity',()=>{
 const poses=['section','entrance','wall'].map(v=>creviceCamera(v,'periodontitis'));assert.equal(new Set(poses.map(p=>JSON.stringify(p))).size,3);
 for(const p of poses)assert.ok([...p.position,...p.target].every(Number.isFinite));
 assert.deepEqual(creviceCamera('wall','periodontitis','td').target,microbialPlacements('periodontitis').find(p=>p.id==='td').position);
});
test('3D examples are independent from immutable saved endpoints and deny measured placement',async()=>{
 const path=new URL('../src/pocket-input.json',import.meta.url),before=await readFile(path);for(const health of ['healthy','periodontitis'])for(const p of buildCrevice(health))p.geometry.dispose();assert.deepEqual(await readFile(path),before);assert.match(CREVICE_CLAIM,/not measured colonization/);assert.match(CREVICE_CLAIM,/sizes exaggerated independently/);
});
test('every authored tissue loft has closed, consistently oriented section boundaries',()=>{
 for(const health of ['healthy','periodontitis'])for(const part of buildCrevice(health)){
  const g=part.geometry;if(part.taxon){g.dispose();continue;}
  const a=g.getAttribute('position'),ix=g.index.array,edges=new Map();let volume=0;
  const key=i=>[a.getX(i),a.getY(i),a.getZ(i)].map(n=>Math.round(n*1e5)).join(',');
  for(let j=0;j<ix.length;j+=3){
   const p=Array.from({length:3},(_,k)=>[a.getX(ix[j+k]),a.getY(ix[j+k]),a.getZ(ix[j+k])]);
   const u=p[1].map((n,k)=>n-p[0][k]),v=p[2].map((n,k)=>n-p[0][k]),cross=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];
   assert.ok(Math.hypot(...cross)>1e-9,health+' '+part.id+' degenerate triangle');
   volume+=p[0][0]*(p[1][1]*p[2][2]-p[1][2]*p[2][1])+p[0][1]*(p[1][2]*p[2][0]-p[1][0]*p[2][2])+p[0][2]*(p[1][0]*p[2][1]-p[1][1]*p[2][0]);
   for(let k=0;k<3;k++){const x=key(ix[j+k]),y=key(ix[j+(k+1)%3]),e=x<y?x+'|'+y:y+'|'+x,count=edges.get(e)||[0,0];count[0]++;count[1]+=x<y?1:-1;edges.set(e,count);}
  }
  assert.ok([...edges.values()].every(v=>v[0]===2&&v[1]===0),health+' '+part.id+' open or reversed seam');assert.ok(volume>0);
  assert.deepEqual(g.groups.map(p=>p.materialIndex),[0,1]);g.dispose();
 }
});
const inside=(pt,poly)=>{let c=false;for(let i=0,j=poly.length-1;i<poly.length;j=i++){const [xi,yi]=poly[i],[xj,yj]=poly[j];if((yi>pt[1])!==(yj>pt[1])&&pt[0]<(xj-xi)*(pt[1]-yi)/(yj-yi)+xi)c=!c;}return c;};
const crosses=(a,b,c,d)=>{const o=(p,q,r)=>Math.sign((q[0]-p[0])*(r[1]-p[1])-(q[1]-p[1])*(r[0]-p[0]));return o(a,b,c)*o(a,b,d)<0&&o(c,d,a)*o(c,d,b)<0;};
test('tissue profiles are simple and do not overlap one another',()=>{
 for(const health of ['healthy','periodontitis']){
  const profiles=creviceProfiles(health),tissues=Object.entries(profiles).filter(([k])=>k!=='cuff');
  for(const [key,p]of Object.entries(profiles)){
   assert.ok(p.length>=8,health+' '+key);
   for(let i=0;i<p.length;i++)for(let j=i+2;j<p.length;j++){if(i===0&&j===p.length-1)continue;assert.ok(!crosses(p[i],p[(i+1)%p.length],p[j],p[(j+1)%p.length]),health+' '+key+' self-intersects');}
  }
  // Dense interior sampling: no point lies inside two tissues.
  let samples=0;for(let x=-6.3;x<2.6;x+=.031)for(let y=-17;y<5.4;y+=.029){const n=tissues.filter(([,p])=>inside([x,y],p)).length;assert.ok(n<=1,health+' overlap at '+x.toFixed(2)+','+y.toFixed(2));samples+=n;}
  assert.ok(samples>20000,health+' profiles cover the section');
 }
});
test('sulcular and junctional epithelia, supracrestal connective tissue and ligament stay distinct',()=>{
 for(const health of ['healthy','periodontitis']){
  const P=creviceProfiles(health),{margin,attachment,junction,crest}=creviceLandmarks(health);
  assert.ok(margin>attachment&&attachment>junction&&junction>crest);
  const near=(key,y,dx)=>inside([toothWall(y)+dx,y],P[key]);
  // Fluid space beside the tooth in the sulcus; junctional epithelium on the tooth below it.
  const mid=(margin+attachment)/2;assert.ok(!Object.entries(P).some(([k,p])=>k!=='cuff'&&inside([toothWall(mid)+creviceGap(mid,health)/2,mid],p)),health+' sulcus is open');
  assert.ok(near('junctional',(attachment+junction)/2,.02),health+' junctional attachment on tooth');
  assert.ok(near('connective',(junction+crest)/2,.05),health+' supracrestal connective tissue meets cementum');
  assert.ok(near('pdl',crest-1,.1)&&near('bone',crest-1,.5),health+' ligament between root and socket');
  assert.ok(near('sulcular',mid,creviceGap(mid,health)+.08),health+' sulcular lining faces the fluid space');
  assert.ok(near('cementum',-2,-.03)&&near('enamel',1,-.05),health+' root and crown coverings');
 }
});
test('section labels point inside the tissue they name, on the opened face',()=>{
 const {axis,radius,half}=CREVICE_SWEEP;
 const profile=([x,y,z])=>[axis+Math.hypot(x-axis,z),y,Math.atan2(z,x-axis)*radius];
 const owner={enamel:'enamel',margin:'gingiva',oral:'gingiva',plaque:'plaque',sulcular:'sulcular',junctional:'junctional',connective:'connective',cementum:'cementum',pdl:'pdl',bone:'bone',dentin:'dentin',pulp:'pulp'};
 for(const health of ['healthy','periodontitis']){
  const P=creviceProfiles(health);
  for(const l of creviceLabels('section',health)){
   const [x,y,z]=profile(l.point);assert.ok(Math.abs(z-half)<1e-6,l.key+' anchor on the section face');
   if(owner[l.key])assert.ok(inside([x,y],P[owner[l.key]]),health+' '+l.key+' anchor outside '+owner[l.key]);
   else assert.ok(!Object.entries(P).some(([k,p])=>k!=='cuff'&&inside([x,y],p)),health+' '+l.key+' anchor should be in the open fluid space');
  }
  assert.ok(creviceLabels('section',health).filter(l=>l.essential).length<=6,'compact layouts keep a short label set');
 }
});
test('illustrative forms rest on the plaque film below the retained margin band',()=>{
 const {axis}=CREVICE_SWEEP;
 for(const health of ['healthy','periodontitis']){
  const P=creviceProfiles(health),cuffBottom=Math.min(...P.cuff.map(p=>p[1]));
  for(const p of microbialPlacements(health)){
   const [x,y,z]=p.position,r=axis+Math.hypot(x-axis,z),gap=r-toothWall(y);
   assert.ok(gap>.05&&gap<.36,health+' '+p.id+' rests on the film, not floating ('+gap.toFixed(3)+')');
   assert.ok(y+(health==='healthy'?.45:.85)*.5<cuffBottom,health+' '+p.id+' visible below the margin band');
   assert.ok(Math.abs(z)<2.4,health+' '+p.id+' inside the swept patch');
  }
 }
});
