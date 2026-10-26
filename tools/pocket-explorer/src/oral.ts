import * as T from 'three';
import {mergeVertices} from 'three/addons/utils/BufferGeometryUtils.js';
import '../tooth-mesh.js';
import '../mouth-geometry.js';
import type {Part} from './geometry';
import type {TissueId} from './anatomy';
declare global {interface Window {ToothMesh:{build:(api:any)=>void};MouthGeometry:Record<string,(api:any)=>void>;}}
// Adapter for existing authored anatomy. Source geometry is preserved; the
// production renderer can later consume reviewed GLB assets through the same IDs.
export function buildMouth():Part[] {
 const parts:Part[]=[];const scale=10;
 for(const upper of [false,true]){
  let material:TissueId='enamel';const buffers=new Map<TissueId,number[]>();
  const map=(p:number[])=>upper?[p[0]*1.04*scale,(2.8-p[1])*scale,p[2]*1.03*scale]:[p[0]*scale,(p[1]-.30)*scale,p[2]*scale];
  const quad=(a:number[],b:number[],c:number[],d:number[])=>{let buf=buffers.get(material);if(!buf){buf=[];buffers.set(material,buf);}for(const p of [a,b,c,a,c,d])buf.push(...map(p));};
  const setMaterial=(m:string)=>{material=({gum:'gingiva',mucosa:upper?'palate':'gingiva',root:'cementum',lip:'gingiva'} as Record<string,TissueId>)[m]||m as TissueId;};
  const center=(a:number)=>[2.55*Math.sin(a),0,3.2*Math.cos(a)];
  const ellipsoid=(c:number[],r:number[])=>{const p=(a:number,b:number)=>[c[0]+r[0]*Math.cos(a)*Math.sin(b),c[1]+r[1]*Math.cos(b),c[2]+r[2]*Math.sin(a)*Math.sin(b)];for(let j=0;j<32;j++)for(let k=0;k<24;k++)quad(p(j/32*Math.PI*2,k/24*Math.PI),p((j+1)/32*Math.PI*2,k/24*Math.PI),p((j+1)/32*Math.PI*2,(k+1)/24*Math.PI),p(j/32*Math.PI*2,(k+1)/24*Math.PI));};
  const api={upper,quad,setMaterial,center,ellipsoid};
  window.MouthGeometry.gingiva(api);
  if(upper)window.MouthGeometry.palate(api);else window.MouthGeometry.tongue(api);
  const angles=[.12,.36,.61,.85,1.07,1.28,1.55];
  for(const side of [-1,1])angles.forEach((a,index)=>window.ToothMesh.build({...api,a:a*side,index,showRoots:true,thirdMode:'absent'}));
  for(const [id,positions]of buffers){const raw=new T.BufferGeometry();raw.setAttribute('position',new T.Float32BufferAttribute(positions,3));const geometry=mergeVertices(raw);raw.dispose();geometry.computeVertexNormals();
   const p=geometry.getAttribute('position'),uv:number[]=[];for(let i=0;i<p.count;i++)uv.push(p.getX(i)/25,p.getZ(i)/25);geometry.setAttribute('uv',new T.Float32BufferAttribute(uv,2));
   if(id==='tongue')geometry.translate(0,9,-1);
   parts.push({id,name:`${upper?'Upper':'Lower'} ${id} context`,geometry,offset:new T.Vector3()});}
 }
 // Concave posterior surround with a genuine opening, rather than a flat sheet.
 const points:number[]=[],uv:number[]=[],indices:number[]=[];const rows=32,n=80;
 for(let i=0;i<=rows;i++)for(let j=0;j<=n;j++){const u=i/rows,a=j/n*Math.PI*2,rx=6+20*u,ry=8+8*u;
  const folds=.7*Math.sin(a*7)*Math.sin(Math.PI*u);
  points.push(rx*Math.cos(a),16+ry*Math.sin(a),-22+22*Math.pow(u,1.7)+folds);uv.push(j/n,u);}
 for(let i=0;i<rows;i++)for(let j=0;j<n;j++){const a=i*(n+1)+j,b=a+n+1;indices.push(a,b,b+1,a,b+1,a+1);}
 const throat=new T.BufferGeometry();throat.setAttribute('position',new T.Float32BufferAttribute(points,3));throat.setAttribute('uv',new T.Float32BufferAttribute(uv,2));throat.setIndex(indices);throat.computeVertexNormals();parts.push({id:'palate',name:'Posterior oral surround',geometry:throat,offset:new T.Vector3()});
 return parts;
}
