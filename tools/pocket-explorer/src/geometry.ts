import * as T from 'three';
import {mergeGeometries} from 'three/addons/utils/BufferGeometryUtils.js';
import type {TissueId} from './anatomy';
export interface Part {id:TissueId;name:string;geometry:T.BufferGeometry;offset:T.Vector3;}
export const CUSPS=[[-2.7,2.5,1.25],[1.8,2.6,1.05],[-2.7,-2.3,1.5],[2,-2.4,1.25],[4.05,.3,.65]];
export function occlusal(x:number,z:number):number {
 let y=5.6;
 for(const [cx,cz,h]of CUSPS)y+=h*Math.exp(-((x-cx)**2+(z-cz)**2)/3.6);
 y-=.26*Math.exp(-z*z*5)*Math.exp(-x*x/22);
 y-=.17*Math.exp(-((x+.18*Math.sin(z*2))**2)*6)*Math.exp(-z*z/15);
 y+=.19*Math.exp(-(((Math.abs(x)-4.25)/.45)**2))*(1-.022*z*z);
 return y;
}
function gridSurface(rows:T.Vector3[][],bottom?:T.Vector3,top?:T.Vector3):T.BufferGeometry {
 const n=rows[0].length,points=rows.flatMap(r=>r.flatMap(v=>v.toArray())),indices:number[]=[];
 for(let r=0;r<rows.length-1;r++)for(let j=0;j<n;j++){const a=r*n+j,b=(r+1)*n+j,c=(r+1)*n+(j+1)%n,d=r*n+(j+1)%n;indices.push(a,b,c,a,c,d);}
 if(bottom){const p=points.length/3;points.push(...bottom.toArray());for(let j=0;j<n;j++)indices.push(p,j,(j+1)%n);}
 if(top){const p=points.length/3;points.push(...top.toArray());const k=(rows.length-1)*n;for(let j=0;j<n;j++)indices.push(k+j,p,k+(j+1)%n);}
 const uv:number[]=[];for(let r=0;r<rows.length;r++)for(let j=0;j<n;j++)uv.push(j/n,r/(rows.length-1));if(bottom)uv.push(.5,0);if(top)uv.push(.5,1);
 const g=new T.BufferGeometry();g.setAttribute('position',new T.Float32BufferAttribute(points,3));g.setAttribute('uv',new T.Float32BufferAttribute(uv,2));g.setIndex(indices);g.computeVertexNormals();return g;
}
export function crownGeometry(scale=1):T.BufferGeometry {
 const n=96,rows:T.Vector3[][]=[],sp=(x:number,p:number)=>Math.sign(x)*Math.abs(x)**p;
 const outline=(t:number)=>[5.45*sp(Math.cos(t),.67)*(1-.07*Math.sin(t)),4.65*sp(Math.sin(t),.76)];
 for(let k=0;k<=30;k++){
  const v=k/30;rows.push(Array.from({length:n},(_,j)=>{
   const [ox,oz]=outline(j/n*Math.PI*2),f=.73+.29*Math.sin(v*Math.PI*.79),x=ox*f,z=oz*f;
   return new T.Vector3(x*scale,(-.45+v*(occlusal(x,z)+.45))*scale,z*scale);
  }));
 }
 const rim=rows.at(-1)!;
 for(let k=1;k<=28;k++){const r=1-k/29;rows.push(rim.map(p=>new T.Vector3(p.x*r,occlusal(p.x*r/scale,p.z*r/scale)*scale,p.z*r)));}
 return gridSurface(rows,new T.Vector3(0,-.45*scale,0),new T.Vector3(0,occlusal(0,0)*scale,0));
}
export function enamelShell():T.BufferGeometry {
 const outer=crownGeometry(),inner=crownGeometry(.84),index=inner.getIndex()!,normals=inner.getAttribute('normal');
 for(let i=0;i<index.count;i+=3){const value=index.getX(i+1);index.setX(i+1,index.getX(i+2));index.setX(i+2,value);}
 for(let i=0;i<normals.count;i++)normals.setXYZ(i,-normals.getX(i),-normals.getY(i),-normals.getZ(i));
 const result=mergeGeometries([outer,inner])!;outer.dispose();inner.dispose();return result;
}
export function rootGeometry(side:number,scale=1):T.BufferGeometry {
 const rows:T.Vector3[][]=[];
 // Divergent tapered roots, flattened mesiodistally, with a shallow longitudinal concavity.
 for(let k=0;k<=50;k++){const v=k/51;rows.push(Array.from({length:56},(_,j)=>{
  const t=j/56*Math.PI*2,taper=Math.pow(1-v,.58),x=side*(2.35+2.3*v)+.24*Math.sin(v*Math.PI),y=-1.6-12.6*v,z=.30*Math.sin(v*Math.PI)*side;
  const flute=1-.12*Math.exp(-(((Math.abs(Math.sin(t))-.18)/.14)**2))*Math.sin(v*Math.PI);
  return new T.Vector3((x+1.58*taper*Math.cos(t)*flute)*scale,y*scale,(z+2.24*taper*Math.sin(t))*scale);
 }));}
 // Reverse vertical sweep winding for the apical direction.
 const g=gridSurface(rows,new T.Vector3(side*2.35*scale,-1.6*scale,0),new T.Vector3(side*4.65*scale,-14.2*scale,0));
 const idx=g.getIndex()!;for(let i=0;i<idx.count;i+=3){const a=idx.getX(i+1);idx.setX(i+1,idx.getX(i+2));idx.setX(i+2,a);}g.computeVertexNormals();return g;
}
export function tube(points:number[][],radius:number):T.BufferGeometry {
 const curve=new T.CatmullRomCurve3(points.map(p=>new T.Vector3(...p as [number,number,number])));
 const g=new T.TubeGeometry(curve,64,radius,12,false);
 // End disks: tangent-normal oriented, closed geometry for sectional contour reconstruction.
 const pos=Array.from(g.getAttribute('position').array),ind=Array.from(g.index!.array);const stride=13;
 for(const [row,t,reverse]of [[0,0,true],[64,1,false]] as const){const center=curve.getPoint(t);const at=pos.length/3;pos.push(...center.toArray());for(let j=0;j<12;j++)ind.push(at,row*stride+(reverse?j+1:j),row*stride+(reverse?j:j+1));}
 g.dispose();const result=new T.BufferGeometry();result.setAttribute('position',new T.Float32BufferAttribute(pos,3));result.setIndex(ind);result.computeVertexNormals();return result;
}
function gumCollar():T.BufferGeometry {
 const rows:T.Vector3[][]=[];
 for(let k=0;k<52;k++){const t=k/52*Math.PI*2;rows.push(Array.from({length:96},(_,j)=>{
  const a=j/96*Math.PI*2;const outer=8.2,inner=4.15,r=inner+(outer-inner)*(.5+.5*Math.cos(t));
  const scallop=.15*Math.pow(Math.abs(Math.cos(a)),8),rise=2.3*(1-(r-inner)/(outer-inner));
  return new T.Vector3(r*Math.sign(Math.cos(a))*Math.abs(Math.cos(a))**.73,-3+rise+.38*Math.sin(t)+scallop,r*.73*Math.sign(Math.sin(a))*Math.abs(Math.sin(a))**.82);
 }));}
 rows.push(rows[0]);return gridSurface(rows);
}
function socketBlock():T.BufferGeometry {
 // Planar polygon extrusion with a socket opening, authored in XY then rotated into XZ.
 const shape=new T.Shape();shape.moveTo(-8,-5.8);shape.lineTo(8,-5.8);shape.quadraticCurveTo(8.8,-5.8,8.8,-4.8);shape.lineTo(8.8,4.8);shape.quadraticCurveTo(8.8,5.8,8,5.8);shape.lineTo(-8,5.8);shape.quadraticCurveTo(-8.8,5.8,-8.8,4.8);shape.lineTo(-8.8,-4.8);shape.quadraticCurveTo(-8.8,-5.8,-8,-5.8);
 const hole=new T.Path();hole.absellipse(0,0,5.5,3.3,0,Math.PI*2,true,0);shape.holes.push(hole);
 const g=new T.ExtrudeGeometry(shape,{depth:11.5,bevelEnabled:true,bevelSegments:3,steps:1,bevelSize:.28,bevelThickness:.28,curveSegments:48});g.rotateX(Math.PI/2);g.translate(0,-3,0);return g;
}
export function buildMolar():Part[] {
 const parts:Part[]=[];const add=(id:TissueId,name:string,g:T.BufferGeometry,offset=[0,0,0])=>parts.push({id,name,geometry:g,offset:new T.Vector3(...offset as [number,number,number])});
 add('enamel','Crown enamel shell',enamelShell(),[5,2,0]);
 add('dentin','Coronal dentin',crownGeometry(.84),[-4,0,0]);
 const trunk=new T.CylinderGeometry(3.9,3.0,2.2,64);trunk.scale(1,1,.83);trunk.translate(0,-1.5,0);add('dentin','Root trunk',trunk,[-4,0,0]);
 for(const side of [-1,1]){
  add('dentin',`${side<0?'Mesial':'Distal'} root dentin`,rootGeometry(side,.94),[-4,0,0]);
  add('cementum','Root cementum',rootGeometry(side),[3,-1,0]);
  const pdl=rootGeometry(side);pdl.scale(1.08,1,1.08);add('pdl','Ligament envelope',pdl,[-3,-2,0]);
 }
 const chamber=new T.SphereGeometry(1,48,32);chamber.scale(2.4,1.4,1.9);chamber.translate(0,2.15,0);add('pulp','Pulp chamber',chamber,[0,0,4]);
 for(const [x,z]of [[-2.2,-1.9],[1.8,-1.9],[-2.2,2],[1.8,2]])add('pulp','Pulp horn',tube([[x,2.1,z*.65],[x,3.35,z],[x*.97,4.5,z]],.3),[0,0,4]);
 const canals=[[-1,0],[-1,.8],[1,0]];
 for(const [side,z]of canals){
  const p=[[side*.8,1.4,z],[side*2.35,-2,z],[side*3.15,-6,z*.8],[side*4.6,-13.95,.08*side]];
  add('pulp','Root canal',tube(p,.34),[0,0,4]);
  add('vessels','Arterial route',tube(p.map(([x,y,z])=>[x+.075,y,z+.055]),.075),[0,0,6]);
  add('vessels','Venous route',tube(p.map(([x,y,z])=>[x-.12,y,z-.05]),.09),[0,0,6]);
  add('nerve','Neural route',tube(p.map(([x,y,z])=>[x,y,z-.18]),.055),[0,0,7]);
 }
 add('gingiva','Gingival collar',gumCollar(),[0,-3,0]);
 add('bone','Alveolar socket block',socketBlock(),[0,-6,0]);
 const septum=new T.Shape();septum.moveTo(-3.55,-14.45);septum.lineTo(3.55,-14.45);septum.quadraticCurveTo(2.5,-7.6,.2,-4.2);septum.quadraticCurveTo(-2.5,-7.6,-3.55,-14.45);
 const septalBone=new T.ExtrudeGeometry(septum,{depth:6.3,bevelEnabled:false,curveSegments:32});septalBone.translate(0,0,-3.15);add('bone','Interradicular bone',septalBone,[0,-6,0]);
 return parts;
}
