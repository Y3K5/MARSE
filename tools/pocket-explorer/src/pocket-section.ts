import * as T from 'three';
import {mergeGeometries} from 'three/addons/utils/BufferGeometryUtils.js';
import {sectionGeometry,type SectionResult} from './section';
// Remove only the anterior right quadrant (x>0 && z>0).
// Intersect closed volumes first, then omit the internal joining face.
export function quarterSection(input:T.BufferGeometry):SectionResult{
 const leftInput=input.clone().rotateY(-Math.PI/2),rightInput=input.clone().rotateY(Math.PI/2);
 const left=sectionGeometry(leftInput,0),right=sectionGeometry(rightInput,0);leftInput.dispose();rightInput.dispose();
 for(const g of [left.surface,left.cap])g.rotateY(Math.PI/2);
 for(const g of [right.surface,right.cap])g.rotateY(-Math.PI/2);
 const rightClosed=mergeGeometries([right.surface,right.cap])!;
 const closed=sectionGeometry(rightClosed,0),rightSkin=sectionGeometry(right.surface,0);
 const rotated=left.cap.clone().rotateY(Math.PI),front=sectionGeometry(rotated,0);front.surface.rotateY(-Math.PI);
 const surface=mergeGeometries([left.surface,rightSkin.surface])!,cap=mergeGeometries([front.surface,closed.cap])!;
 // Object-space projections on each plane keep cut-face surface patterns stable.
 for(const g of [surface,cap]){
  const p=g.getAttribute('position'),n=g.getAttribute('normal'),uv=[];
  for(let i=0;i<p.count;i++)uv.push((Math.abs(n.getX(i))>.8?p.getZ(i):p.getX(i))/12,p.getY(i)/12);
  g.setAttribute('uv',new T.Float32BufferAttribute(uv,2));g.computeBoundingSphere();
 }
 for(const g of [left.surface,left.cap,right.surface,right.cap,rightClosed,closed.surface,closed.cap,rightSkin.surface,rightSkin.cap,rotated,front.surface,front.cap])g.dispose();
 return {surface,cap,loops:left.loops+right.loops+closed.loops,openChains:left.openChains+right.openChains+closed.openChains};
}
