import * as T from 'three';
import {OrbitControls} from '../vendor/OrbitControls.js';
import {sourceNetworkData,fieldLayout,overviewPointGroups} from './exploration-model.js';
import {networkRecords,edgePath,recordStyle} from './network-model.js';

// Independent schematic camera. This never creates a protein correspondence or a source edge.
export class NetworkOverview{
 constructor({canvas,stage,labels,data,colors}){
  this.canvas=canvas;this.stage=stage;this.labels=labels;this.data=sourceNetworkData(data);this.colors=colors;this.dirty=true;
  this.renderer=new T.WebGLRenderer({canvas,antialias:true,alpha:true});this.renderer.setPixelRatio(Math.min(devicePixelRatio,1.5));
  this.camera=new T.PerspectiveCamera(40,1,.1,10000);this.scene=new T.Scene();this.controls=new OrbitControls(this.camera,canvas);this.controls.enableDamping=true;this.controls.enablePan=true;this.controls.minDistance=5;this.controls.maxDistance=5000;
  this.controls.addEventListener('change',()=>this.dirty=true);this.controls.addEventListener('start',()=>this.userCamera=true);
  const taxa=[...data.species.map(s=>s.id)];if([...this.data.proteins,...this.data.contexts].some(p=>p.species==='unknown'))taxa.push('unknown');const centers=Object.fromEntries(taxa.map((id,i)=>[id,[(i-(taxa.length-1)/2)*200,0,(i%2)*12]]));this.layout=fieldLayout(this.data,.7,'cluster',1,centers);this.records=networkRecords(this.data,{scope:'all',candidates:false});this.groupMarks=[];
  this.renderedNodes=0;
  for(const {species,nodes} of overviewPointGroups(this.layout)){
   const sp=data.species.find(s=>s.id===species);
   const positions=nodes.flatMap(n=>n.position),geometry=new T.BufferGeometry();geometry.setAttribute('position',new T.Float32BufferAttribute(positions,3));
   this.scene.add(new T.Points(geometry,new T.PointsMaterial({color:colors[species]||'#b5bdc4',size:5,sizeAttenuation:false})));this.renderedNodes+=nodes.length;
   if(!sp)continue;
   const center=new T.Vector3(...nodes[0].groupCenter),label=document.createElement('span');const fullName=sp.short||sp.name,initials=String(fullName).match(/^([A-Z])[a-z]*\.?\s+([a-z])/i);label.textContent=initials?initials[1].toUpperCase()+initials[2].toLowerCase():fullName;label.title=fullName;label.setAttribute('aria-label',fullName);label.style.color=colors[sp.id];labels.append(label);
   this.groupMarks.push({id:sp.id,center,label,anchor:center.clone().add(new T.Vector3(0,Math.max(...nodes.map(n=>n.position[1]-center.y))+15,0))});
  }
  this.records.forEach((r,i)=>{
   const style=recordStyle(r),curve=new T.CatmullRomCurve3(edgePath(r,this.layout,i).map(p=>new T.Vector3(...p))),geometry=new T.BufferGeometry().setFromPoints(curve.getPoints(36));
   const material=style.pattern==='solid'?new T.LineBasicMaterial({color:style.color,transparent:true,opacity:.8}):new T.LineDashedMaterial({color:style.color,transparent:true,opacity:.65,dashSize:style.pattern==='dot'?.7:2,gapSize:2});
   const line=new T.Line(geometry,material);if(style.pattern!=='solid')line.computeLineDistances();this.scene.add(line);
  });
  this.box=new T.Box3();for(const n of Object.values(this.layout))this.box.expandByPoint(new T.Vector3(...n.position));for(const g of this.groupMarks)this.box.expandByPoint(g.anchor);
  new ResizeObserver(()=>{this.resize();if(!this.userCamera)this.fit();}).observe(stage);this.resize();this.fit();
 }
 resize(){const r=this.stage.getBoundingClientRect();if(!r.width||!r.height)return;this.renderer.setSize(r.width,r.height,false);this.camera.aspect=r.width/r.height;this.camera.updateProjectionMatrix();this.dirty=true;}
 fit(){
  if(this.box.isEmpty())return;this.userCamera=false;
  const center=this.box.getCenter(new T.Vector3()),direction=new T.Vector3(.25,.55,1.7).normalize(),right=new T.Vector3().crossVectors(new T.Vector3(0,1,0),direction).normalize(),up=new T.Vector3().crossVectors(direction,right).normalize(),tan=Math.tan(T.MathUtils.degToRad(this.camera.fov/2));let distance=10;
  for(const x of [this.box.min.x,this.box.max.x])for(const y of [this.box.min.y,this.box.max.y])for(const z of [this.box.min.z,this.box.max.z]){const p=new T.Vector3(x,y,z).sub(center),depth=p.dot(direction);distance=Math.max(distance,depth+(Math.abs(p.dot(up))+16)/tan,depth+(Math.abs(p.dot(right))+16)/(tan*this.camera.aspect));}
  this.controls.enableDamping=false;this.controls.update();this.controls.target.copy(center);this.camera.position.copy(center).addScaledVector(direction,distance*1.08);this.controls.update();this.controls.enableDamping=true;this.dirty=true;
 }
 viewState(){return {position:this.camera.position.toArray(),target:this.controls.target.toArray()};}
 restore(view){this.userCamera=true;this.controls.enableDamping=false;this.controls.update();this.camera.position.fromArray(view.position);this.controls.target.fromArray(view.target);this.controls.update();this.controls.enableDamping=true;this.dirty=true;}
 setActive(species){this.activeSpecies=species;for(const g of this.groupMarks)g.label.classList.toggle('active',g.id===species);this.dirty=true;}
 tick(){
  if(this.stage.closest('[hidden]'))return;this.controls.update();if(!this.dirty)return;
  this.renderer.render(this.scene,this.camera);const r=this.stage.getBoundingClientRect();
  for(const g of this.groupMarks){const p=g.anchor.clone().project(this.camera);g.label.style.left=((p.x+1)/2*r.width)+'px';g.label.style.top=((1-p.y)/2*r.height)+'px';}
  this.canvas.dataset.overviewState=JSON.stringify({species:this.groupMarks.length,nodes:Object.keys(this.layout).length,renderedNodes:this.renderedNodes,records:this.records.length,totalRecords:this.data.relationships.length,activeSpecies:this.activeSpecies,position:this.camera.position.toArray(),target:this.controls.target.toArray()});this.dirty=false;
 }
}
