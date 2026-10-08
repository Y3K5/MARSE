import * as T from 'three';
export interface SectionResult {surface:T.BufferGeometry;cap:T.BufferGeometry;loops:number;openChains:number;}
// Intersect the actual triangles with z <= cut. Reconstruct closed contours and
// triangulate caps with nested holes. No stencil-only or painted cross-section.
export function sectionGeometry(input:T.BufferGeometry,cut:number):SectionResult {
 const g=input.index?input.toNonIndexed():input.clone();
 if(!g.getAttribute('normal'))g.computeVertexNormals();
 const pos=g.getAttribute('position'),norm=g.getAttribute('normal');
 const vertices:number[]=[],normals:number[]=[],segments:Array<[T.Vector3,T.Vector3]>=[];
 const eps=1e-6;
 type V={p:T.Vector3;n:T.Vector3};
 const append=(v:V)=>{vertices.push(...v.p.toArray());normals.push(...v.n.toArray());};
 for(let i=0;i<pos.count;i+=3){
  const tri:V[]=[0,1,2].map(j=>({p:new T.Vector3().fromBufferAttribute(pos,i+j),n:new T.Vector3().fromBufferAttribute(norm,i+j)}));
  const out:V[]=[],hits:T.Vector3[]=[];
  for(let j=0;j<3;j++){
   const a=tri[j],b=tri[(j+1)%3],da=a.p.z-cut,db=b.p.z-cut;
   if(da<=eps)out.push(a);
   if((da>eps&&db<=eps)||(da<=eps&&db>eps)){
    const t=da/(da-db),p=a.p.clone().lerp(b.p,t);p.z=cut;
    const v={p,n:a.n.clone().lerp(b.n,t).normalize()};out.push(v);hits.push(p);
   }
  }
  for(let j=1;j<out.length-1;j++){
   const cross=out[j].p.clone().sub(out[0].p).cross(out[j+1].p.clone().sub(out[0].p));
   if(cross.lengthSq()>1e-16){append(out[0]);append(out[j]);append(out[j+1]);}
  }
  if(hits.length===2&&hits[0].distanceToSquared(hits[1])>1e-14)segments.push([hits[0],hits[1]]);
 }
 const points=new Map<string,T.Vector3>(),edges=new Map<string,Set<string>>();
 const key=(p:T.Vector3)=>`${Math.round(p.x*1e5)},${Math.round(p.y*1e5)}`;
 for(const [a,b]of segments){const ka=key(a),kb=key(b);if(ka===kb)continue;points.set(ka,a);points.set(kb,b);if(!edges.has(ka))edges.set(ka,new Set());if(!edges.has(kb))edges.set(kb,new Set());edges.get(ka)!.add(kb);edges.get(kb)!.add(ka);}
 const loops:T.Vector2[][]=[];let openChains=0;
 while(edges.size){
  const first=edges.keys().next().value as string;let current=first,previous='';const loop:T.Vector2[]=[];let closed=false;
  for(let guard=0;guard<=segments.length+1;guard++){
   const p=points.get(current)!;loop.push(new T.Vector2(p.x,p.y));
   const neighbors=edges.get(current);const next=neighbors&&Array.from(neighbors).find(n=>n!==previous);
   if(!next){edges.delete(current);break;}
   edges.get(current)!.delete(next);edges.get(next)?.delete(current);
   if(edges.get(current)!.size===0)edges.delete(current);
   if(edges.get(next)?.size===0)edges.delete(next);
   previous=current;current=next;
   if(current===first){closed=true;break;}
  }
  if(closed&&loop.length>=3)loops.push(loop);else openChains++;
 }
 const contains=(p:T.Vector2,poly:T.Vector2[])=>{let inside=false;for(let i=0,j=poly.length-1;i<poly.length;j=i++){const a=poly[i],b=poly[j];if((a.y>p.y)!==(b.y>p.y)&&p.x<(b.x-a.x)*(p.y-a.y)/(b.y-a.y)+a.x)inside=!inside;}return inside;};
 const areas=loops.map(l=>Math.abs(T.ShapeUtils.area(l)));
 const parents=loops.map((l,i)=>{let p=-1;for(let j=0;j<loops.length;j++)if(i!==j&&areas[j]>areas[i]&&contains(l[0],loops[j])&&(p<0||areas[j]<areas[p]))p=j;return p;});
 const depth=(i:number):number=>parents[i]<0?0:1+depth(parents[i]);
 const caps:number[]=[];
 loops.forEach((outer,i)=>{
  if(depth(i)%2)return;
  const contour=outer.slice();if(T.ShapeUtils.isClockWise(contour))contour.reverse();
  const holes=loops.filter((_,j)=>parents[j]===i).map(l=>{const h=l.slice();if(!T.ShapeUtils.isClockWise(h))h.reverse();return h;});
  const flat=contour.concat(...holes);
  for(const f of T.ShapeUtils.triangulateShape(contour,holes)){
   const a=flat[f[0]],b=flat[f[1]],c=flat[f[2]];
   const ordered=((b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x))>=0?[a,b,c]:[a,c,b];
   for(const p of ordered)caps.push(p.x,p.y,cut);
  }
 });
 g.dispose();
 const surface=new T.BufferGeometry();surface.setAttribute('position',new T.Float32BufferAttribute(vertices,3));surface.setAttribute('normal',new T.Float32BufferAttribute(normals,3));surface.computeBoundingSphere();
 surface.setAttribute('uv',new T.Float32BufferAttribute(vertices.flatMap((_,i)=>i%3===0?[vertices[i]/10,vertices[i+1]/15]:[]),2));
 const cap=new T.BufferGeometry();cap.setAttribute('position',new T.Float32BufferAttribute(caps,3));cap.computeVertexNormals();cap.computeBoundingSphere();
 cap.setAttribute('uv',new T.Float32BufferAttribute(caps.flatMap((_,i)=>i%3===0?[caps[i]/12,caps[i+1]/12]:[]),2));
 return{surface,cap,loops:loops.length,openChains};
}
