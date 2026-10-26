import * as T from 'three';
import {mergeGeometries,toCreasedNormals} from 'three/addons/utils/BufferGeometryUtils.js';
import {buildPocket,SEGMENT,type PocketPart} from './pocket-geometry';
import {endSection,quarterSection,splitEnds} from './pocket-section';
import type {Health} from './pocket-science';

export interface PreparedPart extends PocketPart {cap?:T.BufferGeometry;}
export interface PreparedView {parts:PreparedPart[];openChains:number;loops:number;triangles:number;}
export const triangleCount=(g:T.BufferGeometry)=>(g.index?.count??g.getAttribute('position').count)/3;

// Cache geometry, not view selections or scientific records. Surfaces belong to
// this cache: removing a Mesh from the scene must not dispose its shared buffers.
export class PocketViewCache {
 private bases=new Map<Health,{parts:PocketPart[];open:number}>();
 private views=new Map<string,PreparedView>();
 builds=0;
 constructor(private decorate?:(part:PreparedPart)=>void){}
 has(health:Health,cutaway:boolean){return this.views.has(health+':'+cutaway);}
 get(health:Health,cutaway:boolean):PreparedView {
  const key=health+':'+cutaway;
  const existing=this.views.get(key);if(existing)return existing;
  let base=this.bases.get(health);
  if(!base){
   let open=0;
   const parts=buildPocket(health).map(p=>{
    const g=toCreasedNormals(p.geometry,T.MathUtils.degToRad(48));p.geometry.dispose();
    const trimmed=endSection(g,SEGMENT.ends[0],SEGMENT.ends[1]);g.dispose();trimmed.cap.dispose();open+=trimmed.openChains;
    return {...p,geometry:trimmed.surface};
   });
   base={parts,open};this.bases.set(health,base);
  }
  let openChains=base.open,loops=0,triangles=0;
  const parts=base.parts.map(p=>{
   let geometry:T.BufferGeometry,cap:T.BufferGeometry;
   if(cutaway){
    const q=quarterSection(p.geometry),s=splitEnds(q.surface,SEGMENT.ends[0],SEGMENT.ends[1]);
    geometry=s.surface;cap=mergeGeometries([q.cap,s.ends])!;openChains+=q.openChains;loops+=q.loops;
    q.surface.dispose();q.cap.dispose();s.ends.dispose();
   }else {const s=splitEnds(p.geometry,SEGMENT.ends[0],SEGMENT.ends[1]);geometry=s.surface;cap=s.ends;}
   const result={...p,geometry,cap};this.decorate?.(result);
   triangles+=triangleCount(geometry)+triangleCount(cap);return result;
  });
  const view={parts,openChains,loops,triangles};this.views.set(key,view);this.builds++;return view;
 }
 dispose(){
  for(const view of this.views.values())for(const p of view.parts){p.geometry.dispose();p.cap?.dispose();}
  for(const base of this.bases.values())for(const p of base.parts)p.geometry.dispose();
  this.views.clear();this.bases.clear();
 }
}
