import * as T from 'three';
import {mergeGeometries} from 'three/addons/utils/BufferGeometryUtils.js';
import {buildMolar,tube,type Part} from './geometry';
import {buildMouth} from './oral';
import {PRESETS,type Health} from './pocket-science';
import organic from './pocket-organic-input.json';
export type PocketId = Part['id']|'plaque'|'lumen'|'epithelium'|'connective'|'supragingival'|'cheek';
export interface PocketPart {id:PocketId;name:string;geometry:T.BufferGeometry;}
export const PROFILES={
 healthy:{lumen:'M4.53,.35 L4.68,.35 Q5.05,-.65 4.85,-1.05 Q4.75,-1.45 4.42,-1.65 L4.40,-1.58 Q4.49,-.60 4.53,.35 Z',
 gingiva:'M4.74,.35 Q5.14,.88 6.05,.28 Q8.8,-.35 9.45,-3.25 L9.75,-7.9 Q7.7,-8.8 5.1,-7.45 L4.5,-3.1 Q5.14,-1.18 4.74,.35 Z'},
 periodontitis:{lumen:'M4.53,-.15 L4.71,-.15 Q5.52,-2.40 5.32,-4.15 Q5.03,-5.68 4.61,-6.20 L4.48,-6.15 Q4.49,-3.9 4.53,-.15 Z',
 gingiva:'M4.80,-.15 Q5.28,.65 6.35,-.18 Q9.1,-1.0 9.66,-4.7 L10.0,-11.1 Q7.5,-11.95 5.14,-10.0 L4.62,-7.0 Q6.08,-3.55 4.80,-.15 Z'}
} as const;
// Parse only our authored M/L/Q/Z profiles, shared with the 2-D fallback.
export function profileShape(profile:string){
 const tokens=profile.match(/[MLQZ]|[-+]?(?:\d*\.)?\d+/g)!;const s=new T.Shape();let i=0;
 while(i<tokens.length){const op=tokens[i++];if(op==='Z'){s.closePath();continue;}const x=+tokens[i++],y=+tokens[i++];if(op==='M')s.moveTo(x,y);else if(op==='L')s.lineTo(x,y);else if(op==='Q')s.quadraticCurveTo(x,y,+tokens[i++],+tokens[i++]);else throw Error('Invalid profile');}return s;
}
function extrude(profile:string,depth:number,bevel=.12){const g=new T.ExtrudeGeometry(profileShape(profile),{depth,bevelEnabled:bevel>0,bevelSize:bevel,bevelThickness:bevel,bevelSegments:3,steps:1,curveSegments:36});g.translate(0,0,-depth/2);return g;}
function shell(outer:T.BufferGeometry,inner:T.BufferGeometry){
 const inside=inner.index?inner.toNonIndexed():inner.clone(),p=inside.getAttribute('position');
 for(let i=0;i<p.count;i+=3)for(let a=0;a<3;a++){const value=p.array[(i+1)*3+a];p.array[(i+1)*3+a]=p.array[(i+2)*3+a];p.array[(i+2)*3+a]=value;}
 const outside=outer.index?outer.toNonIndexed():outer.clone();inside.deleteAttribute('uv');outside.deleteAttribute('uv');inside.deleteAttribute('normal');outside.deleteAttribute('normal');
 const g=mergeGeometries([outside,inside])!;g.computeVertexNormals();outside.dispose();inside.dispose();return g;
}
export function buildScaffoldPocket(health:Health):PocketPart[]{
 const base=buildMolar(),parts:PocketPart[]=[];
 for(const p of base){
  if(['gingiva','bone','pdl'].includes(p.id))continue;
  if(p.id==='cementum'){const inner=p.geometry.clone();inner.scale(.94,.94,.94);parts.push({id:p.id,name:p.name,geometry:shell(p.geometry,inner)});inner.dispose();}
  else parts.push({id:p.id,name:p.name,geometry:p.geometry.clone()});
 }
 for(const p of base.filter(p=>p.id==='cementum')){
  const outer=p.geometry.clone();outer.scale(1.07,1.008,1.07);parts.push({id:'pdl',name:'Periodontal ligament shell',geometry:shell(outer,p.geometry)});outer.dispose();
 }
 for(const p of base)p.geometry.dispose();
 const v=PRESETS[health];
 parts.push({id:'connective',name:'Gingival connective tissue',geometry:extrude(PROFILES[health].gingiva,7.0)});
 // Surface cover, separated from the cut-face connective tissue.
 const cover=extrude(PROFILES[health].gingiva,7.30);parts.push({id:'gingiva',name:'Gingival surface',geometry:cover});
 parts.push({id:'lumen',name:'Pocket fluid space',geometry:extrude(PROFILES[health].lumen,3.8,0)});
 const e=`M4.73,${v.margin} Q5.72,-2.0 ${health==='healthy'?'4.99,-1.10':'5.47,-4.18'} Q5.10,${v.attachment-.25} 4.62,${v.attachment-.75} L4.47,${v.attachment-.72} Q5.02,${v.attachment-.18} ${health==='healthy'?'4.83,-1.05':'5.31,-4.12'} Q5.41,-2.0 4.59,${v.margin} Z`;
 // Healthy lining has a separate short profile; no disease shape is reused.
 const healthyE='M4.69,.35 L4.84,.35 Q5.20,-.65 4.98,-1.10 Q4.84,-1.60 4.53,-2.4 L4.40,-2.36 Q4.72,-1.45 4.84,-1.02 Q5.00,-.60 4.69,.35 Z';
 parts.push({id:'epithelium',name:'Sulcular and junctional epithelium',geometry:extrude(health==='healthy'?healthyE:e,4.0,0)});
 parts.push({id:'plaque',name:'Tooth-attached plaque',geometry:extrude(`M4.45,${v.margin-.20} L4.52,${v.margin-.20} L4.59,${v.attachment+.15} L4.49,${v.attachment+.12} Z`,4.15,0)});
 const leftG='M-4.74,.35 Q-5.4,.80 -6.4,.10 Q-9.1,-.8 -9.5,-3.6 L-9.6,-8.8 Q-6.3,-9.1 -4.9,-7.8 L-4.4,-2.5 Q-4.8,-1.1 -4.74,.35 Z';
 parts.push({id:'gingiva',name:'Opposing gingival context',geometry:extrude(leftG,7.3)});
 const rightB=`M10.2,${v.crest} Q7.2,${v.crest+.75} 5.02,${v.crest} Q5.48,-11.0 5.02,-14.5 Q3.3,-15.5 2.8,-14.7 L2.0,-16.25 Q5.8,-17.4 9.0,-16.1 Q11.2,-13.8 10.2,${v.crest} Z`;
 const leftB='M-10.2,-3.5 Q-7.2,-2.8 -5.03,-3.5 Q-5.50,-11.0 -5.03,-14.5 Q-3.3,-15.5 -2.8,-14.7 L-2.0,-16.25 Q-5.8,-17.4 -9.0,-16.1 Q-11.2,-13.8 -10.2,-3.5 Z';
 const septum='M-2.75,-15.6 L2.75,-15.6 Q2.7,-10.2 .05,-4.1 Q-2.7,-10.2 -2.75,-15.6 Z';
 for(const [name,profile]of [['Selected-site alveolar bone',rightB],['Opposing alveolar bone',leftB],['Interradicular bone',septum]])parts.push({id:'bone',name,geometry:extrude(profile,6.6,.10)});
 parts.push({id:'vessels',name:'Gingival vascular route',geometry:tube([[8.3,-10,-.8],[7.7,-7,-.8],[7.1,-4,-.6],[6.2,-2,-.5]],.085)});
 return parts;
}
export function buildOralContext():PocketPart[]{
 const all=buildMouth();for(const p of all.filter(p=>p.id==='cementum'))p.geometry.dispose();
 const parts:PocketPart[]=all.filter(p=>p.id!=='cementum').map(p=>({id:p.id,name:p.name,geometry:p.geometry}));
 for(const p of parts)if(p.id==='tongue')p.geometry.translate(0,-5,0);
 // Author a contextual bony arch below each gingival arch; no bone measurements.
 for(const upper of [false,true]){
  const arc=new T.CatmullRomCurve3(Array.from({length:81},(_,i)=>{const a=-1.75+i/80*3.5;return new T.Vector3(25.5*Math.sin(a),upper?39:-10,32*Math.cos(a));}));
  const g=new T.TubeGeometry(arc,112,3.7,16,false);parts.push({id:'bone',name:upper?'Maxillary arch context':'Mandibular arch context',geometry:g});
 }
 return parts;
}
// Three-tooth context segment recorded by the Blender build (FDI 35, 36, 37).
export const SEGMENT=(()=>{const s=organic.segment;if(s.ends.length!==2||!(s.ends[0]<s.ends[1]))throw Error('Incompatible segment ends');return {...s,ends:[s.ends[0],s.ends[1]] as [number,number]};})();
const smooth=(t:number)=>{const u=Math.min(1,Math.max(0,t));return u*u*(3-2*u);};
// Authored oral-side context for FDI 36. A lower first molar faces the cheek
// (buccal, +z) and the tongue (lingual, -z); lips border the front teeth.
// Each sheet starts under the gum's basal edge as alveolar mucosa, folds at the
// vestibule or floor of the mouth, and fades toward the specimen ends.
function sideSheet(inner:Array<[number,number]>,thickness:number,side:1|-1,bend:number,x0:number,x1:number,alpha:(x:number,y:number)=>number){
 // Offset the inner profile along its normal, away from the tooth.
 const outer=inner.map(([z,y],i)=>{const a=inner[Math.max(i-1,0)],b=inner[Math.min(i+1,inner.length-1)],tz=b[0]-a[0],ty=b[1]-a[1],l=Math.hypot(tz,ty)||1;let nz=ty/l,ny=-tz/l;if(nz*side<0){nz=-nz;ny=-ny;}return [z+thickness*nz,y+thickness*ny] as [number,number];});
 const s=new T.Shape();
 s.moveTo(...inner[0]);for(const p of inner.slice(1))s.lineTo(...p);for(const p of outer.reverse())s.lineTo(...p);s.closePath();
 const g=new T.ExtrudeGeometry(s,{depth:x1-x0,steps:36,bevelEnabled:true,bevelSize:.12,bevelThickness:.25,bevelSegments:3,curveSegments:24});
 const p=g.getAttribute('position'),colors:number[]=[];for(let i=0;i<p.count;i++){const u=p.getX(i),v=p.getY(i),x=p.getZ(i)+x0;p.setXYZ(i,x,v,u+bend*x*x);colors.push(1,1,1,alpha(x,v));}
 g.setAttribute('color',new T.Float32BufferAttribute(colors,4));g.deleteAttribute('uv');g.computeVertexNormals();return g;
}
export function buildOralSides():PocketPart[]{
 const [x0,x1]=SEGMENT.ends,bend=.006;
 const curve=(pts:Array<[number,number]>)=>new T.SplineCurve(pts.map(([z,y])=>new T.Vector2(z,y))).getPoints(48).map(v=>[v.x,v.y] as [number,number]);
 const ends=(x:number,a:number,b:number)=>smooth(Math.min(x-a,b-x)/2.2);
 // Cheek side: mucosa leaves the gum edge, folds at the buccal vestibule and
 // rises as the cheek lining. Drawn over the mesial part so the pocket stays clear.
 const cheekEnd=-.6;
 const cheek=sideSheet(curve([[5.75,-7.2],[5.95,-8.9],[6.8,-10.3],[8.1,-9.7],[8.9,-6],[9.3,0],[9.5,5.2]]),.7,1,bend,x0,cheekEnd,(x,y)=>.34*ends(x,x0,cheekEnd)*(1-.55*smooth((y+1)/6)));
 // Tongue side: floor of the mouth runs from the lingual gum edge to the tongue.
 const floor=sideSheet(curve([[-5.75,-7.2],[-5.95,-8.9],[-6.8,-10.3],[-8.2,-9.6],[-9.3,-6.5],[-10.3,-.9]]),.7,-1,bend,x0,x1,(x)=>.40*ends(x,x0,x1));
 const tongue=new T.SphereGeometry(1,64,40);tongue.scale(12.5,3.4,4.6);
 const tp=tongue.getAttribute('position'),tc:number[]=[];for(let i=0;i<tp.count;i++){const x=tp.getX(i)+.45;tp.setXYZ(i,x,tp.getY(i)+1.6-.006*x*x,tp.getZ(i)-12.3+.02*x*x);tc.push(1,1,1,.55*ends(x,x0,x1));}
 tongue.setAttribute('color',new T.Float32BufferAttribute(tc,4));tongue.computeVertexNormals();
 return [{id:'cheek',name:'Buccal mucosa, vestibule and cheek',geometry:cheek},{id:'tongue',name:'Floor of the mouth (tongue side)',geometry:floor},{id:'tongue',name:'Tongue lateral border',geometry:tongue}];
}
export const ANCHORS:Partial<Record<PocketId,[number,number,number]>>={supragingival:[-1.6,1.1,4.1],cheek:[-6,-2,9.4],tongue:[-6,3,-7.8],enamel:[-2.55,5.9,2.05],dentin:[-1.6,1.5,0],pulp:[.7,2.2,0],cementum:[4.25,-8.5,0],pdl:[4.7,-9,0],gingiva:[-5.12,-.25,2.7],connective:[7.5,-5,0],bone:[7.4,-8.4,0],plaque:[4.5,-3,0],lumen:[4.8,-3.5,0],epithelium:[4.9,-6.1,0]};

export function buildPocket(health:Health):PocketPart[]{
 if(organic.format!=='marse.pocket-organic-packed/2'||organic.calibration!==null)throw Error('Incompatible authored geometry');
 return organic.cases[health].map(index=>{
  const part=organic.pool[index],decode=(s:string)=>Uint8Array.from(atob(s),c=>c.charCodeAt(0));
  const positions=new Float32Array(decode(part.positions).buffer),indices=new Uint32Array(decode(part.indices).buffer);
  const g=new T.BufferGeometry();g.setAttribute('position',new T.BufferAttribute(positions,3));g.setIndex(new T.BufferAttribute(indices,1));g.computeVertexNormals();
  const p=g.getAttribute('position'),uv:number[]=[];for(let i=0;i<p.count;i++)uv.push(p.getX(i)/18,p.getY(i)/26);g.setAttribute('uv',new T.Float32BufferAttribute(uv,2));
  return {id:part.id as PocketId,name:part.name,geometry:g};
 });
}
