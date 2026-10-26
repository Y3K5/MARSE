import * as T from 'three';
import atlas from './pocket-atlas-input.json';
import type {ContextPart} from './pocket-context';
import type {PocketId} from './pocket-geometry';

export const ATLAS_SOURCE=atlas.source;
export const ATLAS_CREDIT=atlas.attribution;
// Source coordinates are mm. This rigid basis change preserves the common
// coordinate system. The jaw opening below is authored, not measured TMJ motion.
export const ATLAS_BASIS=new T.Matrix4().set(1,0,0,0, 0,0,1,-1460, 0,-1,0,-140, 0,0,0,1);
const pivot=new T.Vector3(0,40,-45);
export const ATLAS_OPEN=new T.Matrix4().makeTranslation(...pivot.toArray()).multiply(new T.Matrix4().makeRotationX(.45)).multiply(new T.Matrix4().makeTranslation(...pivot.clone().negate().toArray()));
export interface AtlasPart extends ContextPart {atlasIdentity:{element:string;concept:string;name:string;fdi:number|null;sourceSHA256:string};}
function floatBuffer(encoded:string){const bytes=Uint8Array.from(atob(encoded),c=>c.charCodeAt(0));return new Float32Array(bytes.buffer);}
function indexBuffer(encoded:string){const bytes=Uint8Array.from(atob(encoded),c=>c.charCodeAt(0));return new Uint32Array(bytes.buffer);}
export function atlasContext(scope:'mouth'|'face'):AtlasPart[]{
 return atlas.parts.filter(p=>(scope==='face'||p.scope==='mouth')&&p.id!=='vessels'&&p.id!=='connective').map(p=>{
  const g=new T.BufferGeometry();g.setAttribute('position',new T.BufferAttribute(floatBuffer(p.positions),3));g.setIndex(new T.BufferAttribute(indexBuffer(p.indices),1));
  g.applyMatrix4(ATLAS_BASIS);if(p.jaw==='lower')g.applyMatrix4(ATLAS_OPEN);g.computeVertexNormals();
  // UVs supply artistic material detail only. They are not histology maps.
  const a=g.getAttribute('position'),uv=[];
  for(let i=0;i<a.count;i++)uv.push(a.getX(i)/60+.5,(a.getY(i)+a.getZ(i)*.2)/70+.5);
  g.setAttribute('uv',new T.Float32BufferAttribute(uv,2));
  return {id:p.id as PocketId,name:p.fdi?`FDI ${p.fdi} · atlas whole tooth`:p.name,geometry:g,
   atlasIdentity:{element:p.element,concept:p.concept,name:p.name,fdi:p.fdi,sourceSHA256:p.source_sha256}};
 });
}
const selected=atlasContext('mouth').find(p=>p.atlasIdentity.fdi===36)!;
selected.geometry.computeBoundingBox();
// Marker is a display anchor on the identified atlas tooth. It does not register
// the independent three-tooth pocket specimen onto this atlas.
export const ATLAS_FDI36=selected.geometry.boundingBox!.getCenter(new T.Vector3()).add(new T.Vector3(0,4,0));
selected.geometry.dispose();
