import * as T from 'three';
import {OrbitControls} from '../vendor/OrbitControls.js';
import {makeFold} from './fold-geometry.js';
import {regionResidues} from './interaction-model.js';
import {regionOverlay} from './interaction-geometry.js';
import {siteResidueEntries,siteLabelCandidates} from './site-residue-model.js';
export class InteractionLens{
 constructor({canvas,stage,labels=null,onResidue}){
  Object.assign(this,{canvas,stage,labelHost:labels,onResidue,dirty:true,siteRows:[],labelItems:[]});
  this.renderer=new T.WebGLRenderer({canvas,antialias:true,alpha:true});this.renderer.setPixelRatio(Math.min(devicePixelRatio,1.5));this.renderer.outputColorSpace=T.SRGBColorSpace;
  this.scene=new T.Scene();this.scene.add(new T.HemisphereLight('#e4f3ef','#233542',2.4));const light=new T.DirectionalLight('#ffe3bb',3);light.position.set(30,60,80);this.scene.add(light);
  this.camera=new T.PerspectiveCamera(40,1,.05,10000);this.controls=new OrbitControls(this.camera,canvas);this.controls.enableDamping=true;this.controls.minDistance=.2;this.controls.maxDistance=3000;this.controls.addEventListener('change',()=>this.dirty=true);
  this.observer=new ResizeObserver(()=>this.resize());this.observer.observe(stage);let down;
  canvas.addEventListener('pointerdown',e=>down=[e.clientX,e.clientY]);
  canvas.addEventListener('pointerup',e=>{
   if(!down||Math.hypot(e.clientX-down[0],e.clientY-down[1])>5||!this.fold)return;
   const rect=canvas.getBoundingClientRect(),ray=new T.Raycaster();ray.setFromCamera(new T.Vector2((e.clientX-rect.left)/rect.width*2-1,-(e.clientY-rect.top)/rect.height*2+1),this.camera);
   // Gold overlays are intentionally visible over faded context, so they get the same picking priority.
   const overlays=this.fold.children.filter(c=>c.userData.regionIndices),context=this.fold.children.filter(c=>c.userData.residueInstances),hit=ray.intersectObjects(overlays,false)[0]||ray.intersectObjects(context,false)[0];
   if(hit?.instanceId!==undefined)this.selectResidue(hit.object.userData.regionIndices?.[hit.instanceId]??hit.instanceId);
  });
  canvas.addEventListener('keydown',e=>{
   if(!['+','=','-','Home','ArrowLeft','ArrowRight','ArrowUp','ArrowDown'].includes(e.key))return;e.preventDefault();if(e.key==='Home'){this.fit(this.focusBox);return;}
   const offset=this.camera.position.clone().sub(this.controls.target),s=new T.Spherical().setFromVector3(offset);
   if(e.key==='ArrowLeft')s.theta-=.14;if(e.key==='ArrowRight')s.theta+=.14;if(e.key==='ArrowUp')s.phi=Math.max(.05,s.phi-.14);if(e.key==='ArrowDown')s.phi=Math.min(Math.PI-.05,s.phi+.14);if(['+','='].includes(e.key))s.radius*=.8;if(e.key==='-')s.radius*=1.25;
   this.camera.position.copy(this.controls.target).add(new T.Vector3().setFromSpherical(s));this.controls.update();this.dirty=true;
  });
 }
 clear(){
  this.dirty=true;this.siteRows=[];this.selectedIndex=null;this.selectionMarker=null;this.labelItems=[];this.labelHost?.replaceChildren();
  if(this.fold){this.scene.remove(this.fold);this.fold.traverse(o=>{o.geometry?.dispose();o.material?.dispose();});this.fold=null;}
 }
 show(data,ledger,location,region){
  this.clear();this.location=location;this.region=region;
  const model=region?.modelProteinId||location.modelProteinId||location.endpointId,shape=data.shape_traces[model];
  if(!shape){this.canvas.dataset.lensState=JSON.stringify({endpoint:location.endpointId,shape:false,exactContactResolved:false,regions:0});return;}
  this.fold=makeFold(shape,ledger.entries[model],{actual:true,quality:'high',colorMode:'path',beads:true});this.scene.add(this.fold);let regionBox=null;
  if(region){
   this.siteRows=siteResidueEntries(data,ledger,location,region);
   for(const child of this.fold.children){child.material.vertexColors=false;child.material.color.set('#66818b');if(child.isInstancedMesh&&child.instanceColor){for(let i=0;i<child.count;i++)child.setColorAt(i,new T.Color('#66818b'));child.instanceColor.needsUpdate=true;}child.material.transparent=true;child.material.opacity=child.isInstancedMesh?.12:.22;child.material.depthWrite=false;child.material.needsUpdate=true;}
   const rows=regionResidues(data,ledger,location.recordId,location.endpointId,region.id),overlay=regionOverlay(shape,rows,{actual:true,annotationKind:region.category});overlay.material.depthTest=false;overlay.material.depthWrite=false;overlay.renderOrder=20;this.fold.add(overlay);regionBox=new T.Box3().setFromObject(overlay).expandByScalar(4);
  }
  this.fold.updateMatrixWorld(true);this.wholeBox=new T.Box3().setFromObject(this.fold);this.focusBox=regionBox||this.wholeBox;this.buildLabels();this.resize();this.fit(this.focusBox);this.dirty=true;
 }
 selectResidue(index,{focus=false}={}){
  const row=this.fold?.userData.residues[index];if(!row)return;
  this.selectedIndex=index;const entry=this.siteRows.find(r=>r.index===index)||null;
  if(this.selectionMarker){this.selectionMarker.removeFromParent();this.selectionMarker.geometry.dispose();this.selectionMarker.material.dispose();this.selectionMarker=null;}
  if(entry){
   const marker=new T.Mesh(new T.SphereGeometry(1.45,14,10),new T.MeshBasicMaterial({color:'#fff2ad',wireframe:true,transparent:true,opacity:.8,depthTest:false,depthWrite:false}));marker.position.set(...row.point.slice(1).map((v,k)=>v-this.fold.userData.sourceCenter[k]));marker.renderOrder=25;this.fold.add(marker);this.selectionMarker=marker;
  }
  this.buildLabels();this.onResidue(row,{entry,region:this.region});if(focus)this.focusSelected();this.dirty=true;
 }
 focusSelected(){
  if(!this.fold||this.selectedIndex==null)return;const row=this.fold.userData.residues[this.selectedIndex],point=new T.Vector3(...row.point.slice(1).map((v,k)=>v-this.fold.userData.sourceCenter[k]));this.fit(new T.Box3(point.clone().addScalar(-4),point.clone().addScalar(4)));
 }
 buildLabels(){
  if(!this.labelHost)return;this.labelHost.replaceChildren();this.labelItems=[];
  this.labelLines=document.createElementNS('http://www.w3.org/2000/svg','svg');this.labelLines.setAttribute('aria-hidden','true');this.labelHost.append(this.labelLines);
  for(const row of siteLabelCandidates(this.siteRows,this.selectedIndex)){
   const button=document.createElement('button');button.className='site-label';button.dataset.siteLabel=String(row.index);button.textContent=row.label;button.title=`${row.fullName} (${row.name}) · paper ${row.paper} → canonical ${row.canonical} → author ${row.authorLabel}`;button.setAttribute('aria-label',button.title);button.setAttribute('aria-pressed',String(row.index===this.selectedIndex));button.onclick=()=>this.selectResidue(row.index);this.labelHost.append(button);this.labelItems.push({row,button});
  }
 }
 updateLabels(){
  if(!this.labelHost||!this.fold)return;const rect=this.stage.getBoundingClientRect(),placed=[],lines=[];this.labelLines?.setAttribute('viewBox',`0 0 ${rect.width} ${rect.height}`);
  for(const {row,button}of this.labelItems){
   const p=new T.Vector3(...row.point.slice(1).map((v,k)=>v-this.fold.userData.sourceCenter[k])).applyMatrix4(this.fold.matrixWorld).project(this.camera),x=(p.x+1)/2*rect.width,y=(1-p.y)/2*rect.height;
   let chosen=null;const w=button.offsetWidth||90,h=button.offsetHeight||25;
   if(p.z>-1&&p.z<1&&x>3&&x<rect.width-3&&y>3&&y<rect.height-3){
    for(const [dx,dy]of [[14,-h-10],[14,12],[-w-14,-h-10],[-w-14,12],[-w/2,-h-34],[-w/2,34]]){
     const box={x:Math.max(4,Math.min(rect.width-w-4,x+dx)),y:Math.max(4,Math.min(rect.height-h-4,y+dy)),w,h};
     if(!placed.some(b=>box.x<b.x+b.w+4&&box.x+box.w>b.x-4&&box.y<b.y+b.h+4&&box.y+box.h>b.y-4)){chosen=box;break;}
    }
   }
   button.hidden=!chosen;if(!chosen)continue;placed.push(chosen);button.style.left=chosen.x+'px';button.style.top=chosen.y+'px';
   const endX=Math.max(chosen.x,Math.min(chosen.x+chosen.w,x)),endY=Math.max(chosen.y,Math.min(chosen.y+chosen.h,y));lines.push(`<line x1="${x}" y1="${y}" x2="${endX}" y2="${endY}"/>`);
  }
  if(this.labelLines)this.labelLines.innerHTML=lines.join('');
 }
 zoomBy(factor){const delta=this.camera.position.clone().sub(this.controls.target).multiplyScalar(factor);this.camera.position.copy(this.controls.target).add(delta);this.controls.update();this.dirty=true;}
 fit(box){
  if(!box||box.isEmpty())return;const center=box.getCenter(new T.Vector3()),direction=new T.Vector3(.6,.35,1).normalize(),right=new T.Vector3().crossVectors(new T.Vector3(0,1,0),direction).normalize(),up=new T.Vector3().crossVectors(direction,right).normalize(),tan=Math.tan(T.MathUtils.degToRad(this.camera.fov/2));let distance=5;
  for(const x of [box.min.x,box.max.x])for(const y of [box.min.y,box.max.y])for(const z of [box.min.z,box.max.z]){const p=new T.Vector3(x,y,z).sub(center),depth=p.dot(direction);distance=Math.max(distance,depth+(Math.abs(p.dot(up))+2)/tan,depth+(Math.abs(p.dot(right))+2)/(tan*this.camera.aspect));}
  this.controls.enableDamping=false;this.controls.update();this.camera.position.copy(center).addScaledVector(direction,distance*1.08);this.controls.target.copy(center);this.controls.update();this.controls.enableDamping=true;this.dirty=true;
 }
 resize(){const {width,height}=this.stage.getBoundingClientRect();if(width<1||height<1)return;this.renderer.setSize(width,height,false);this.camera.aspect=width/height;this.camera.updateProjectionMatrix();this.dirty=true;}
 tick(){
  if(this.stage.closest('[hidden]'))return;this.controls.update();if(!this.dirty)return;this.renderer.render(this.scene,this.camera);this.updateLabels();
  const entry=this.siteRows.find(r=>r.index===this.selectedIndex);
  this.canvas.dataset.lensState=JSON.stringify({endpoint:this.location?.endpointId,record:this.location?.recordId,shape:!!this.fold,region:this.region?{...this.region}:null,exactContactResolved:false,highlightedResidues:this.fold?.children.find(o=>o.userData.regionIndices)?.count||0,siteResidues:this.siteRows.length,selectedResidue:entry?{name:entry.name,chain:entry.chain,author:entry.author,insertion:entry.insertion,paper:entry.paper,canonical:entry.canonical,annotationId:entry.annotationId}:null,visibleLabels:this.labelItems.filter(l=>!l.button.hidden).map(l=>({index:l.row.index,text:l.row.label})),position:this.camera.position.toArray(),target:this.controls.target.toArray(),residueCount:this.fold?.userData.coordinateCount||0});this.dirty=false;
 }
}
