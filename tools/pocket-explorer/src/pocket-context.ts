import * as T from 'three';
import {ParametricGeometry} from 'three/addons/geometries/ParametricGeometry.js';
import {tube} from './geometry';
import type {PocketPart} from './pocket-geometry';
import face from './pocket-face-input.json';

// Authored shared oral coordinates, not atlas registration. +x is the subject's
// left, +z anterior, +y superior. The detail frame has +x distal and +z buccal.
export const FDI36_POSE={angle:1.28,position:[25.5*Math.sin(1.28),0,32*Math.cos(1.28)] as [number,number,number]};
export const SEGMENT_TO_MOUTH=new T.Matrix4().makeRotationY(Math.atan2(32*Math.sin(FDI36_POSE.angle),25.5*Math.cos(FDI36_POSE.angle))).setPosition(...FDI36_POSE.position);
export interface ContextPart extends PocketPart {appearance?:'skin'|'lip'|'eye'|'iris'|'vein'|'artery'|'nerve';}
const ellipsoid=(x:number,y:number,z:number,rx:number,ry:number,rz:number)=>{const g=new T.SphereGeometry(1,64,40);g.scale(rx,ry,rz);g.translate(x,y,z);return g;};
function arc(radius:number,y:number,depth:number,thickness:number){return tube(Array.from({length:41},(_,i)=>{const a=-1.73+i/40*3.46;return [radius*Math.sin(a),y,depth*Math.cos(a)];}),thickness);}
// Context crown: named tooth classes, anatomical surface cues, not measured teeth.
function crown(index:number,upper:boolean){
 const width=[2.8,3.0,3.65,3.75,4.15,5.15,4.85][index],depth=[2.3,2.4,3.1,3.5,3.65,4.5,4.45][index];
 const cusps=index<2?[[-.6,0,.15],[.6,0,.15]]:index===2?[[0,0,2]]:index<5?[[0,1.65,1.5],[0,-1.65,.85]]:index===6?[[-2,2,1.2],[2,2,1],[-2,-2,1.45],[2,-2,1.35]]:[[-2,2,1.2],[1,2,1.15],[3,.8,.75],[-2,-2,1.5],[1.8,-2,1.4]];
 const g=new ParametricGeometry((u,v,p)=>{
  const a=u*Math.PI*2,r=Math.sin(v*Math.PI),cx=Math.sign(Math.cos(a))*Math.abs(Math.cos(a))**.72,cz=Math.sign(Math.sin(a))*Math.abs(Math.sin(a))**.76;
  const x=width*r*cx,z=depth*r*cz,cap=Math.max(0,-Math.cos(v*Math.PI));
  let y=3.2-3.2*Math.cos(v*Math.PI);
  if(v>.5){for(const [qx,qz,h] of cusps)y+=h*Math.exp(-((x-qx)**2+(z-qz)**2)/2.5)*cap;y-=index>=3?.25*Math.exp(-z*z*2)*cap:0;}
  p.set(x,upper?-y:y,z);
 },64,42);
 // Upper crowns reverse the vertical axis. Orient each closed crown outwards.
 const ids=g.index!,p=g.getAttribute('position');let volume=0;
 for(let i=0;i<ids.count;i+=3){const a=new T.Vector3().fromBufferAttribute(p,ids.getX(i)),b=new T.Vector3().fromBufferAttribute(p,ids.getX(i+1)),c=new T.Vector3().fromBufferAttribute(p,ids.getX(i+2));volume+=a.dot(b.cross(c));}
 if(volume<0)for(let i=0;i<ids.count;i+=3){const b=ids.getX(i+1);ids.setX(i+1,ids.getX(i+2));ids.setX(i+2,b);}
 g.computeVertexNormals();return g;
}
export function oralContext():ContextPart[]{
 const parts:ContextPart[]=[];
 for(const upper of [false,true]){
  const cej=upper?31:0,angles=[.12,.36,.61,.85,1.07,1.28,1.55];
  for(const side of [-1,1])for(let index=0;index<7;index++){
   // The lower left first molar is inserted from the actual detail asset later.
   if(!upper&&side===1&&index===5)continue;
   const a=angles[index]*side,position=[25.5*Math.sin(a)*(upper?1.04:1),cej,32*Math.cos(a)*(upper?1.03:1)] as [number,number,number];
   const g=crown(index,upper);g.rotateY(Math.atan2(32*Math.sin(a),25.5*Math.cos(a)));g.translate(...position);
   const tooth=upper?(side===1?21:11)+index:(side===1?31:41)+index;
   parts.push({id:'enamel',name:`FDI ${tooth} · authored ${index<2?'incisor':index===2?'canine':index<5?'premolar':'molar'} crown`,geometry:g});
  }
  parts.push({id:'gingiva',name:upper?'Upper gingival arch':'Lower gingival arch',geometry:arc(25.5,cej+(upper?2.0:-2),32,3.3)});
  parts.push({id:'bone',name:upper?'Maxillary alveolar arch':'Mandibular alveolar arch',geometry:arc(25.5,cej+(upper?7.1:-8.2),32,4.7)});
 }
 // Curved palate with an anterior roof and posterior soft-tissue continuation.
 const palate=new ParametricGeometry((u,v,p)=>{const x=(u-.5)*42,z=-13+v*43;p.set(x,34.5-4*(x/21)**2-2*Math.sin(v*Math.PI),z);},64,56);
 parts.push({id:'palate',name:'Hard palate · authored vault',geometry:palate});
 parts.push({id:'palate',name:'Soft palate · posterior continuation',geometry:ellipsoid(0,29,-13,19,3.4,7)});
 parts.push({id:'palate',name:'Uvula · schematic',geometry:ellipsoid(0,24,-17,1.5,4.5,1.6)});
 const tongue=ellipsoid(0,1,9,20,6,22),pos=tongue.getAttribute('position');
 for(let i=0;i<pos.count;i++){const x=pos.getX(i),y=pos.getY(i),z=pos.getZ(i);pos.setY(i,y-.38*Math.exp(-x*x/2)*Math.max(0,(y+2)/7)+.045*Math.sin(x*3)*Math.sin(z*2));}tongue.computeVertexNormals();
 parts.push({id:'tongue',name:'Tongue · authored open-mouth pose',geometry:tongue});
 parts.push({id:'tongue',name:'Floor of mouth',geometry:ellipsoid(0,-7,7,22,2.6,23)});
 for(const side of [-1,1]){
  parts.push({id:'cheek',name:side===1?'Left buccal mucosa and vestibule':'Right buccal mucosa and vestibule',geometry:ellipsoid(side*34,14,9,3.0,16,23)});
  const ramus=tube([[side*24,-8,-8],[side*34,-5,-16],[side*35,17,-23],[side*33,35,-25]],4.1);
  parts.push({id:'bone',name:'Mandibular ramus · authored context',geometry:ramus});
  parts.push({id:'bone',name:'Mandibular condyle · schematic',geometry:ellipsoid(side*33,38,-25,6,3,3)});
 }
 return parts;
}
export function facialContext():ContextPart[]{
 const skin=new T.BufferGeometry();skin.setAttribute('position',new T.Float32BufferAttribute(face.positions,3));skin.setIndex(face.indices);skin.computeVertexNormals();
 const p=skin.getAttribute('position'),uv:number[]=[];
 for(let i=0;i<p.count;i++)uv.push(p.getX(i)/160+.5,p.getY(i)/240);
 skin.setAttribute('uv',new T.Float32BufferAttribute(uv,2));
 const parts:ContextPart[]=[{id:'cheek',name:'MakeHuman CC0 facial context · authored registration',geometry:skin,appearance:'skin'}];
 for(const side of [-1,1]){
  // Registered stock graphical landmarks, not measured orbital anatomy.
  parts.push({id:'cheek',name:'Eye · graphical facial context',geometry:ellipsoid(side*24.616,90.378,25.624,8.2,8.2,8.2),appearance:'eye'});
  parts.push({id:'cheek',name:'Iris · graphical facial context',geometry:ellipsoid(side*24.616,90.378,33.6,3.1,3.1,.4),appearance:'iris'});
 }
 return parts;
}
// These coordinates represent the named pathway, not measured axon trajectories.
export function neuralContext(whole:boolean):ContextPart[]{
 const parts:ContextPart[]=[];
 const add=(name:string,pts:number[][],r:number,appearance:ContextPart['appearance'])=>parts.push({id:appearance==='nerve'?'nerve':'vessels',name,geometry:tube(pts,r),appearance});
 for(const [x,apex,z]of [[-1.6,-14.45,-.12],[-1.6,-14.45,.12],[3.9,-13.6,0]]){
  for(const [shift,r,style]of [[.11,.065,'nerve'],[-.10,.075,'artery'],[-.26,.09,'vein']] as const){
   const side=x<0?-1:1,pts=[[x+shift,-17,z],[x+shift,apex-.2,z],[side*2.85+shift,-9,z],[side*2.0+shift,-3,z],[side*.65+shift,.8,z]];
   const g=tube(pts,r);if(whole)g.applyMatrix4(SEGMENT_TO_MOUTH);
   parts.push({id:style==='nerve'?'nerve':'vessels',name:`Pulpal ${style} · enlarged schematic diameter`,geometry:g,appearance:style});
  }
 }
 if(whole){
  const c=FDI36_POSE.position;
  add('Inferior alveolar nerve → V3 → sensory root (schematic)',[[c[0],-17,c[2]],[27,-15,-4],[31,-6,-17],[32,18,-22],[25,49,-24],[17,66,-22],[9,75,-17],[4,72,-24]],.5,'nerve');
  add('Mental branch · lower lip and chin',[[27,-15,-4],[23,-14,19],[15,-18,33],[8,-16,40]],.30,'nerve');
  add('Inferior alveolar arterial context',[[c[0]-.4,-16.2,c[2]],[27,-14.2,-4],[30,0,-18]],.35,'artery');
  add('Inferior alveolar venous context',[[c[0]+.45,-17.5,c[2]],[28,-15.5,-4],[32,0,-18]],.40,'vein');
 }else{
  add('Inferior alveolar dental nerve · local schematic',[[ -12,-17,0],[0,-17,0],[12,-17,0]],.16,'nerve');
  add('Arterial route · local schematic',[[-12,-16.5,-.1],[0,-16.5,-.1],[12,-16.5,-.1]],.12,'artery');
  add('Venous route · local schematic',[[-12,-17.5,-.3],[0,-17.5,-.3],[12,-17.5,-.3]],.15,'vein');
 }
 return parts;
}
