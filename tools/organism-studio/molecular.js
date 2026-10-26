// Local atom-coordinate rendering; no network services and no molecular inference.
const colors=['#a9d5cd','#d9a97a','#a7bde7','#d7b6d4','#c4d68f','#e7bdaf'];
export class MolecularStudio{
 constructor(element,onSelect){this.element=element;this.onSelect=onSelect;this.viewer=null;this.model=null;this.chains=[];this.current=null;this.representation='cartoon';this.count=999;this.explode=0;this.surface=null;this.token=0;this.atomCount=0;this.showContacts=false;this.contactCount=0;}
 async load(record,frame='source'){
  const token=++this.token;this.current=record;this.frame=frame;
  let url=`structures/${record.id}.${record.format||'pdb'}`,format=record.format||'pdb';this.placement=null;this.fit=null;
  if(frame==='organism'){url=`integration/placements/${record.id}.pdb`;format='pdb';this.placement=await fetch(`integration/placements/${record.id}.json`).then(r=>{if(!r.ok)throw Error('No Organism Orient placement for this model');return r.json()});}
  if(frame==='memorient'){url=`integration/memorient/${record.id}/oriented.pdb`;format='pdb';this.fit=await fetch(`integration/memorient/${record.id}/result.json`).then(r=>{if(!r.ok)throw Error('No MemOrient fit for this model');return r.json()});}
  const raw=await fetch(url).then(r=>{if(!r.ok)throw Error('Coordinate file unavailable');return r.text()});if(token!==this.token)return false;
  if(!this.viewer)this.viewer=$3Dmol.createViewer(this.element,{backgroundColor:'#0e1b24',antialias:true,disableFog:false});
  this.viewer.clear();this.model=this.viewer.addModel(raw,format,{keepH:true,doAssembly:false});
  const allowed=new Set(Object.keys(record.chains));this.model.removeAtoms(this.model.selectedAtoms({}).filter(a=>!a.hetflag&&!allowed.has(a.chain)));const atoms=this.model.selectedAtoms({});this.atomCount=atoms.length;
  for(const a of atoms){a._original=[a.x,a.y,a.z];}
  this.chains=[...new Set(atoms.filter(a=>!a.hetflag&&a.atom==='CA').map(a=>a.chain))];this.count=this.chains.length;this.explode=0;
  this.apply();this.viewer.zoomTo({hetflag:false});if(frame!=='source')this.viewer.rotate(65,'x');this.viewer.zoom(1.3);this.viewer.resize();this.viewer.render();return true;
 }
 apply(){if(!this.model)return;
  const atoms=this.model.selectedAtoms({}),active=this.chains.slice(0,this.count),center=[0,0,0];let n=0;const centers={};
  for(const a of atoms){a.x=a._original[0];a.y=a._original[1];a.z=a._original[2];if(a.atom==='CA'){n++;for(let k=0;k<3;k++)center[k]+=a._original[k];const c=centers[a.chain]||(centers[a.chain]=[0,0,0,0]);for(let k=0;k<3;k++)c[k]+=a._original[k];c[3]++;}}
  for(let k=0;k<3;k++)center[k]/=Math.max(1,n);
  for(const a of atoms){if(this.explode&&centers[a.chain]){const c=centers[a.chain],v=c.slice(0,3).map((x,k)=>x/c[3]-center[k]),len=Math.hypot(...v)||1;a.x+=v[0]/len*this.explode;a.y+=v[1]/len*this.explode;a.z+=v[2]/len*this.explode;}}
  this.viewer.removeAllSurfaces();this.viewer.removeAllShapes();this.highlightShape=null;this.viewer.removeAllLabels();this.model.setStyle({},{});
  for(const [i,c]of active.entries()){const col=colors[(this.current.chainEntities?.[c]?Object.keys(this.current.entities).indexOf(this.current.chainEntities[c]):this.chains.indexOf(c))%colors.length];let style;
   if(this.representation==='stick')style={stick:{radius:.16,colorscheme:'Jmol'},cartoon:{color:col,opacity:.13}};
   else if(this.representation==='backbone')style={cartoon:{style:'trace',color:col,thickness:.35}};
   else style={cartoon:{color:col,thickness:.4,opacity:this.representation==='surface'?.25:1}};
   this.model.setStyle({chain:c},style);}
  this.model.setStyle({hetflag:true,resn:['HOH','WAT']},{});
  this.viewer.addStyle({hetflag:true,invert:false,not:{resn:['HOH','WAT']}},{stick:{radius:.13,colorscheme:'Jmol'}});
  if(this.representation==='surface'){const selected=this.model.selectedAtoms({chain:active});if(selected.length>25000)throw Error('Choose fewer chains before generating a molecular surface (25,000 atom limit).');this.viewer.addSurface($3Dmol.SurfaceType.MS,{opacity:.5,color:'#b4d5c9'},{chain:active});}
  this.model.setClickable({},false);this.model.setClickable({atom:'CA',chain:active},true,(atom)=>{this.select(atom.chain,String(atom.resi)+(atom.icode||''));this.onSelect(atom);});
  let membrane=null;
  if(this.placement)membrane={half:this.placement.surface_spec.thickness/2,size:this.placement.surface_spec.size};
  if(this.fit)membrane={half:this.fit.summary.half_thickness,size:[140,140]};
  if(membrane){const [w,h]=membrane.size;for(const z of [-membrane.half,membrane.half]){this.viewer.addBox({center:{x:0,y:0,z},dimensions:{w,h,d:.8},color:z>0?'#d1b186':'#7eafbf',opacity:.3});}this.viewer.addLabel(this.fit?'MemOrient hydrophobic boundaries':'Organism Orient geometric interfaces',{position:{x:-w/2,y:-h/2,z:membrane.half},fontSize:11,backgroundOpacity:.65,fontColor:'#e5efed'});}
  const visible=this.model.selectedAtoms({chain:active,atom:'CA'});if(visible.length){const xs=visible.map(a=>a.x),ys=visible.map(a=>a.y),zs=visible.map(a=>a.z);const x=Math.min(...xs),y=Math.min(...ys)-10,z=Math.min(...zs);const length=20;this.viewer.addLabel(`${this.current.id} · ${this.frame} · ${this.explode?'exploded display':'source-linked geometry'}`,{position:{x,y:Math.max(...ys)+8,z},fontSize:12,backgroundOpacity:.7,fontColor:'#dbe7e7'});
   this.viewer.addLine({start:{x,y,z},end:{x:x+length,y,z},color:'#ffffff',linewidth:2});this.viewer.addLabel('2 nm · 20 Å',{position:{x:x+length/2,y:y-3,z},fontSize:11,backgroundOpacity:0,fontColor:'#dbe7e7'});}
  this.contactCount=0;
  if(this.showContacts){const cas=this.model.selectedAtoms({atom:'CA',hetflag:false,chain:active}),grid=new Map(),pairs=new Map(),key=p=>p.map(x=>Math.floor(x/8));for(const a of cas){const k=key(a._original).join(',');if(!grid.has(k))grid.set(k,[]);grid.get(k).push(a);}for(const a of cas){const g=key(a._original);for(let x=-1;x<=1;x++)for(let y=-1;y<=1;y++)for(let z=-1;z<=1;z++){for(const b of grid.get([g[0]+x,g[1]+y,g[2]+z].join(','))||[]){if(a.chain===b.chain)continue;const d=a._original.reduce((s,v,k)=>s+(v-b._original[k])**2,0);if(d>64)continue;const id=[a.chain,b.chain].sort().join(':');if(!pairs.has(id)||pairs.get(id).d>d)pairs.set(id,{a,b,d});}}}this.contactCount=pairs.size;for(const {a,b} of pairs.values()){this.viewer.addLine({start:{x:a.x,y:a.y,z:a.z},end:{x:b.x,y:b.y,z:b.z},color:'#f1dd9d',linewidth:2,dashed:true});}this.viewer.addLabel(`${pairs.size} inter-chain Cα neighbors ≤8 Å · not chemical bonds`,{position:{x:Math.min(...cas.map(a=>a.x)),y:Math.max(...cas.map(a=>a.y))+15,z:Math.min(...cas.map(a=>a.z))},fontSize:11,backgroundOpacity:.6,fontColor:'#f1dd9d'});}
  this.viewer.render();
 }
 select(chain,resid){if(!this.viewer)return;const atom=this.model.selectedAtoms({chain,resi:parseInt(resid,10),atom:'CA'})[0];if(!atom)return;if(this.highlightShape)this.viewer.removeShape(this.highlightShape);this.highlightShape=this.viewer.addSphere({center:{x:atom.x,y:atom.y,z:atom.z},radius:1.25,color:'#ffffff',opacity:.9});this.viewer.render();}

 resize(){this.viewer?.resize();}
 reset(){if(!this.viewer)return;this.viewer.zoomTo({hetflag:false});this.viewer.zoom(1.3);this.viewer.render();}
 screenshot(){if(!this.viewer)return;const old=this.viewer.getView(),w=this.element.clientWidth,h=this.element.clientHeight;this.viewer.setWidth(2400);this.viewer.setHeight(1800);this.viewer.setView(old);this.viewer.render();const uri=this.viewer.pngURI();this.viewer.setWidth(w);this.viewer.setHeight(h);this.viewer.setView(old);this.viewer.resize();this.viewer.render();const a=document.createElement('a');a.href=uri;a.download=`${this.current.id}-${this.frame}-${this.representation}.png`;a.click();}
 info(){return{atoms:this.atomCount,chains:this.chains.length,visibleChains:Math.min(this.count,this.chains.length),representation:this.representation,frame:this.frame,explodedAngstrom:this.explode,proximityLinks:this.contactCount,engineFit:this.fit?.summary,placement:this.placement?.scientific_state};}
}
