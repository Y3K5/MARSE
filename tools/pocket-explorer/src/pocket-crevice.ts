import * as T from 'three';
import {mergeGeometries,mergeVertices} from 'three/addons/utils/BufferGeometryUtils.js';
import {TAXA,PRESETS,type Health} from './pocket-science';
import {MORPHOLOGY} from './pocket-morphology';
import type {PocketId} from './pocket-geometry';
export type CreviceView='section'|'entrance'|'wall';
export type CreviceFinish='enamel'|'hard'|'bone'|'pulp'|'ligament'|'soft'|'lining'|'film'|'cell';
export interface CrevicePart {id:PocketId;name:string;geometry:T.BufferGeometry;color:string;taxon?:string;source?:string;detail?:string;finish?:CreviceFinish;views?:CreviceView[]}
export const CREVICE_CLAIM='Magnified authored crevice; widths and cell sizes exaggerated independently. Example placements are not measured colonization, contacts or density-derived cells.';
type P=[number,number];
// Every tissue is a closed (x, y) profile swept around one tapering root axis:
// x runs from the tooth outward, y runs crown-up, and z is arc length along the
// root surface at the cementoenamel junction (y = 0). Sweeping the same profile
// makes the gingiva wrap the root as a cuff. These are authored display
// coordinates, not a reconstructed root or a physical-scale basis.
const AXIS=-6.4,CORE=AXIS+.2,RADIUS=-1.65-AXIS,HALF=3.2,STEPS=30,LATTICE=.14;
export const CREVICE_DEPTH=HALF*2;
/** Sweep constants shared with the renderer's edge fade. */
export const CREVICE_SWEEP={axis:AXIS,radius:RADIUS,half:HALF};
const ss=(x:number,a:number,b:number)=>T.MathUtils.smoothstep(x,a,b);
// Dentin surface: root taper and a slight apical curvature below the
// cementoenamel junction; a gentle cervical bulge beneath the enamel above it.
function dentinWall(y:number){return y<=0?-1.65+.045*y-.0011*y*y:-1.65+.14*Math.sin(Math.min(y/2.6,1)*Math.PI/2);}
// Cementum thickens apically; enamel thins to a feather edge at the junction.
function cementum(y:number){return .035+.11*(1-ss(y,-3,0))+.004*Math.max(0,-y);}
function enamelOuter(y:number){const t=Math.min(Math.max(y,0)/2.6,1);return -1.615+.43*Math.sin(t*Math.PI/2)**.85;}
const canal=(y:number)=>AXIS+1.3+.035*y;
/** Outer tooth surface: enamel above the cementoenamel junction, cementum below it. */
export const toothWall=(y:number)=>y<0?dentinWall(y)+cementum(y):enamelOuter(y);
export function crevicePoint(x:number,y:number,z:number):[number,number,number]{
 const r=Math.max(x-AXIS,.05),a=z/RADIUS;return [AXIS+r*Math.cos(a),y,r*Math.sin(a)];
}
/** Outward unit normal of the swept root surface at arc position z. */
export const creviceRadial=(z:number):[number,number,number]=>[Math.cos(z/RADIUS),0,Math.sin(z/RADIUS)];
export function creviceLandmarks(health:Health){
 const p=PRESETS[health],margin=p.margin*1.6,attachment=p.attachment*1.6,crest=p.crest*1.6;
 // The junctional length (half the attachment-to-crest distance) is an authored
 // proportion that keeps supracrestal connective tissue visible, not a measurement.
 return {margin,attachment,crest,junction:attachment-.5*(attachment-crest),ligament:crest-.2,bottom:crest-3.6};
}
export function creviceGap(y:number,health:Health){
 const {margin,attachment}=creviceLandmarks(health),t=T.MathUtils.clamp((margin-y)/(margin-attachment),0,1);
 // A narrow entrance under the margin, an exaggerated interior, closing to the
 // junctional attachment at the base.
 return .16*(1-t)**1.5+(health==='healthy'?.78:1.45)*Math.sin(Math.PI*t**1.1)**.9;
}
// Points on a global y lattice so neighbouring tissues that follow the same
// curve receive identical vertices along their shared boundary.
function lattice(a:number,b:number){
 const out=[a],up=b>a;let k=up?Math.floor(a/LATTICE)+1:Math.ceil(a/LATTICE)-1;
 for(;up?k*LATTICE<b-1e-6:k*LATTICE>b+1e-6;k+=up?1:-1)if(Math.abs(k*LATTICE-a)>1e-6)out.push(k*LATTICE);
 out.push(b);return out;
}
const along=(f:(y:number)=>number,a:number,b:number):P[]=>lattice(a,b).map(y=>[f(y),y]);
function spline(points:P[],n:number):P[]{return new T.CatmullRomCurve3(points.map(([x,y])=>new T.Vector3(x,y,0)),false,'centripetal').getPoints(n).map(v=>[v.x,v.y]);}
function chain(...parts:P[][]):P[]{
 const out:P[]=[];for(const part of parts)for(const p of part){const q=out.at(-1);if(!q||Math.hypot(q[0]-p[0],q[1]-p[1])>1e-7)out.push(p);}
 const f=out[0],l=out.at(-1)!;if(Math.hypot(f[0]-l[0],f[1]-l[1])<1e-7)out.pop();return out;
}
// Shared landmarks and boundary curves for one preset.
function frame(health:Health){
 const L=creviceLandmarks(health),{margin:M,attachment:A,junction:J,crest:C,ligament:PT,bottom:B}=L,ts=toothWall;
 const S=(y:number)=>ts(y)+creviceGap(y,health);
 const sulcular=(y:number)=>.22+.5*Math.exp(-(y-A)/.8);
 const junctional=(y:number)=>.06+.66*((y-J)/(A-J))**1.3;
 // Epithelium-to-connective-tissue boundary of the sulcular and junctional epithelia.
 const E=(y:number)=>y>=A?S(y)+sulcular(y):ts(y)+junctional(y);
 const sM=S(M),eM=E(M);
 // Outer (oral) gingival surface with a faint free-gingival groove.
 const O=(y:number)=>{const u=M+.15-y;return sM+.62+2*(1-Math.exp(-u/1.4))+.045*u-.07*Math.exp(-(((u-1.9)/.45)**2));};
 const thM=O(M)-eM;
 // Oral epithelium inner boundary; its wave is an artistic rete-ridge cue.
 const G=(y:number)=>{const v=M-y;return O(y)-(.3+(thM-.3)*Math.exp(-v/.6)+.06*Math.sin(v*7)*(1-Math.exp(-v/.8)));};
 const pdl=.3,boneOuter=(y:number)=>ts(y)+pdl+1.1+.38*(C-.6-y);
 const plaqueTop=M+.3,plaqueBottom=A+.35,cuffDepth=Math.min(.9,.16*(M-A));
 const film=(y:number)=>.015+.085*Math.max(0,Math.min(1,(plaqueTop-y)/.35,(y-plaqueBottom)/.6))**.6+.01*Math.sin(y*4.3)**2;
 const arc=spline([[sM,M],[sM+.03,M+.16],[sM+.14,M+.29],[sM+.32,M+.34],[sM+.5,M+.27],[O(M+.15),M+.15]],22);
 const crestArc=spline([[ts(PT)+pdl,PT],[ts(C)+pdl+.08,C-.07],[ts(C)+pdl+.32,C],[ts(C)+pdl+.66,C-.1],[ts(C)+pdl+.9,C-.42],[boneOuter(C-.6),C-.6]],20);
 return {L,M,A,J,C,PT,B,ts,S,E,O,G,sM,eM,pdl,boneOuter,plaqueTop,plaqueBottom,cuffDepth,film,arc,crestArc,junctional};
}
const e26=enamelOuter(2.6),d26=dentinWall(2.6);
const ENAMEL_TOP:P[]=[[e26,2.6],[e26-.08,3.4],[e26-.5,4.35],[e26-1.3,5],[e26-2.2,5.15],[e26-3,4.75],[AXIS+1,4.35],[CORE,4.25]];
const DEJ_TOP:P[]=[[d26,2.6],[d26-.12,3.3],[d26-.6,3.95],[d26-1.4,4.3],[d26-2.2,4.25],[AXIS+1.1,3.75],[CORE,3.6]];
const CHAMBER:P[]=[[canal(0),0],[AXIS+1.6,.8],[AXIS+1.65,1.7],[AXIS+1.2,2.25],[AXIS+.6,2.3],[CORE,2.2]];
/** Closed profiles for each authored tissue; exported for geometry checks. */
export function creviceProfiles(health:Health):Record<string,P[]>{
 const f=frame(health),{M,A,J,PT,B,ts,S,E,O,G,pdl,boneOuter,plaqueTop,plaqueBottom,film,arc,crestArc}=f,C=f.C;
 const dejTop=spline(DEJ_TOP,40),chamber=spline(CHAMBER,28);
 return {
  enamel:chain([[dentinWall(0),0]],along(enamelOuter,0,2.6),spline(ENAMEL_TOP,48),[...dejTop].reverse(),along(dentinWall,2.6,0)),
  dentin:chain(along(dentinWall,B,2.6),dejTop,[[CORE,3.6]],[...chamber].reverse(),along(canal,0,B)),
  pulp:chain([[CORE,B]],along(canal,B,0),chamber),
  cementum:chain(along(dentinWall,0,B),along(ts,B,0)),
  sulcular:chain(along(S,M,A),along(E,A,M)),
  junctional:chain(along(ts,A,J),along(E,J,A)),
  gingiva:chain(arc,along(O,M+.15,B),along(G,B,M)),
  connective:chain(along(E,M,A),along(E,A,J),[[ts(J),J]],along(ts,J,PT),crestArc,along(boneOuter,C-.6,B),along(G,B,M)),
  pdl:chain(along(ts,PT,B),along(y=>ts(y)+pdl,B,PT)),
  bone:chain(crestArc,along(boneOuter,C-.6,B),along(y=>ts(y)+pdl,B,PT)),
  plaque:chain(along(ts,plaqueTop,plaqueBottom),along(y=>ts(y)+film(y),plaqueBottom,plaqueTop)),
  // Retained marginal band for the biofilm-surface view only.
  cuff:chain(along(S,M-f.cuffDepth,M),arc,along(O,M+.15,M-f.cuffDepth))
 };
}
function loft(source:P[],relief?:(x:number,y:number,z:number)=>number){
 const ring=source.map(([x,y])=>new T.Vector2(x,y));
 // Collinear samples along straight closures carry no silhouette information
 // and can make the cap triangulator emit zero-area triangles.
 for(let i=ring.length-1;i>=0;i--){const a=ring[(i+ring.length-1)%ring.length],b=ring[i],c=ring[(i+1)%ring.length];if(Math.abs((b.x-a.x)*(c.y-b.y)-(b.y-a.y)*(c.x-b.x))<1e-10)ring.splice(i,1);}
 if(T.ShapeUtils.isClockWise(ring))ring.reverse();
 const n=ring.length;
 // Profile corners sharper than 50° get separate vertex columns so shading
 // does not smear across them; smooth samples share one column.
 const sharp=ring.map((b,i)=>{const a=ring[(i+n-1)%n],c=ring[(i+1)%n],u=b.clone().sub(a),v=c.clone().sub(b);return u.dot(v)/(u.length()*v.length())<Math.cos(50*Math.PI/180);});
 const into:number[]=[],out:number[]=[],column:number[]=[];
 for(let i=0;i<n;i++){into[i]=column.push(i)-1;out[i]=sharp[i]?column.push(i)-1:into[i];}
 const length:number[]=[0];for(let i=1;i<n;i++)length[i]=length[i-1]+ring[i].distanceTo(ring[i-1]);
 const positions:number[]=[],uv:number[]=[],indices:number[]=[];
 // Optional outward relief vanishes at both ends so the swept surface still
 // meets its flat section caps exactly.
 for(let j=0;j<=STEPS;j++){const z=-HALF+2*HALF*j/STEPS;for(const i of column){const p=ring[i],lift=relief&&j>0&&j<STEPS?relief(p.x,p.y,z):0;positions.push(...crevicePoint(p.x+lift,p.y,z));uv.push(length[i]/3,z/3);}}
 const width=column.length;
 for(let j=0;j<STEPS;j++)for(let i=0;i<n;i++){
  const a=j*width+out[i],b=j*width+into[(i+1)%n],c=a+width,d=b+width;indices.push(a,b,c,b,d,c);
 }
 const surfaceCount=indices.length,triangles=T.ShapeUtils.triangulateShape(ring,[]);
 // Separate cap vertices keep a clean normal break and a dedicated section
 // material. Their positions exactly match the swept end rings.
 for(const end of [0,STEPS]){
  const offset=positions.length/3,z=-HALF+2*HALF*end/STEPS;
  for(const p of ring){positions.push(...crevicePoint(p.x,p.y,z));uv.push((p.x+7)/3.2,(p.y+17)/3.2);}
  for(const [a,b,c]of triangles)indices.push(...(end===0?[offset+a,offset+c,offset+b]:[offset+a,offset+b,offset+c]));
 }
 // Orient the entire closed volume consistently outward.
 let volume=0;for(let i=0;i<indices.length;i+=3){const a=indices[i]*3,b=indices[i+1]*3,c=indices[i+2]*3;volume+=positions[a]*(positions[b+1]*positions[c+2]-positions[b+2]*positions[c+1])+positions[a+1]*(positions[b+2]*positions[c]-positions[b]*positions[c+2])+positions[a+2]*(positions[b]*positions[c+1]-positions[b+1]*positions[c]);}
 if(volume<0)for(let i=0;i<indices.length;i+=3)[indices[i+1],indices[i+2]]=[indices[i+2],indices[i+1]];
 const g=new T.BufferGeometry();g.setAttribute('position',new T.Float32BufferAttribute(positions,3));g.setAttribute('uv',new T.Float32BufferAttribute(uv,2));g.setIndex(indices);g.addGroup(0,surfaceCount,0);g.addGroup(surfaceCount,indices.length-surfaceCount,1);g.computeVertexNormals();g.computeBoundingBox();g.userData.closedLoft=true;return g;
}
function merged(parts:T.BufferGeometry[]){const flat=parts.map(g=>g.index?g.toNonIndexed():g),out=mergeGeometries(flat,false)!;for(const g of new Set([...parts,...flat]))g.dispose();out.computeVertexNormals();return out;}
function sphere(x:number,y:number,z:number,r:number){const g=new T.SphereGeometry(r,32,24);g.translate(x,y,z);return g;}
function path(points:number[][],radius:number){const curve=new T.CatmullRomCurve3(points.map(p=>new T.Vector3(...p as [number,number,number])));return merged([new T.TubeGeometry(curve,96,radius,12,false),sphere(...points[0] as [number,number,number],radius),sphere(...points.at(-1)! as [number,number,number],radius)]);}
// Qualitative body morphology only. Counts, sizes, spacing and surface relief are
// authored display choices, independent of saved fields and spatial observations.
export function microbialObject(id:string):T.BufferGeometry{
 const m=MORPHOLOGY[id];if(!m)throw Error('Unknown microbial object');let g:T.BufferGeometry;
 if(m.form==='cocci')g=id==='ss'?merged(Array.from({length:5},(_,i)=>sphere(.04*Math.sin(i),-.64+i*.32,.025*Math.cos(i),.185))):id==='vp'?merged([sphere(0,-.22,0,.245),sphere(.035,.23,.04,.245)]):merged([sphere(0,-.23,0,.255),sphere(.025,.26,.06,.23)]);
 else if(m.form==='branching')g=merged([path([[0,-.85,0],[0,-.3,.03],[0,.25,-.02],[0,.85,.06]],.105),path([[0,-.15,0],[.01,.18,.27],[0,.5,.51]],.09),path([[0,.12,0],[-.015,.45,-.25],[0,.65,-.47]],.083)]);
 else if(m.form==='spirochete')g=path(Array.from({length:65},(_,i)=>{const t=i/64;return [.07*Math.sin(t*19+.6)+.035*Math.sin(t*39),-1.2+2.4*t,.17*Math.sin(t*26)+.05*Math.sin(t*11)];}),.042);
 else if(m.form==='fusiform'||m.form==='long-rod'){
  const length=m.form==='fusiform'?2.05:1.75,radius=m.form==='fusiform'?.18:.17;
  const pts=Array.from({length:33},(_,i)=>{const t=i/32;return new T.Vector2(radius*Math.sin(Math.PI*t)**(m.form==='fusiform'?1.2:.4),-length/2+length*t);});g=new T.LatheGeometry(pts,40);
 }else g=new T.CapsuleGeometry(id==='pg'?.235:.18,id==='pg'?.43:.7,12,32);
 // Barely visible authored surface variation, not envelope ultrastructure.
 if(!g.index){g.deleteAttribute('normal');g.deleteAttribute('uv');const smooth=mergeVertices(g,1e-5);g.dispose();g=smooth;}
 const a=g.getAttribute('position');for(let i=0;i<a.count;i++){const x=a.getX(i),y=a.getY(i),z=a.getZ(i),r=1+.009*Math.sin(y*21)*Math.cos(z*17);a.setXYZ(i,x*r,y,z*r);}g.computeVertexNormals();g.computeBoundingBox();return g;
}
// Each body lies flat in its local y-z plane; this is its half-thickness along
// local x, used to rest it on the plaque film rather than float in the fluid.
const REST:Record<string,number>={ss:.19,ac:.12,vp:.25,fn:.18,pi:.18,pg:.235,tf:.17,td:.12,pm:.26};
const cellScale=(health:Health)=>health==='healthy'?.45:.85;
export function microbialPlacements(health:Health){
 const f=frame(health),span=f.M-f.A-f.cuffDepth;
 // All identities use the same three display rows below the retained margin
 // band. No succession or depth niche is assigned. Rows are a gallery
 // arrangement, not positions from a specimen.
 return TAXA.map((t,i)=>{const row=Math.floor(i/3),y=f.M-f.cuffDepth-span*(.16+row*.29),z=-1.9+(i%3)*1.9;return {id:t.id,position:crevicePoint(f.ts(y)+f.film(y)+REST[t.id]*cellScale(health)-.015,y,z),z};});
}
const DETAILS={
 enamel:['Enamel','Mineralized crown surface; it thins to a feather edge at the cementoenamel junction.','Coronal context drawn as a smooth dome. Thickness and contour are authored, not a measured enamel map.'],
 dentin:['Root dentin','Mineralized tissue forming the bulk of the root around the pulp.','The root tapers and curves slightly toward the apex in this drawing. Dentinal tubules are not resolved.'],
 pulp:['Pulp','Living connective tissue in the pulp chamber and root canal.','A tapering canal and chamber roof as authored context; canal geometry varies between real teeth.'],
 cementum:['Cementum','Thin mineralized layer covering the root; ligament and supracrestal fibres insert into it.','Drawn thickening apically from a feather edge at the junction. Thickness is exaggerated and authored.'],
 plaque:['Attached plaque','Tooth-attached microbial biofilm, above and below the gingival margin, distinct from the fluid space.','A thin film on the tooth surface ending short of the junctional attachment. Thickness is illustrative and no species positions are assigned to it.'],
 lumen:['Sulcus · fluid space','The space between the tooth surface and the sulcular epithelium, open at the gingival margin.','Width is exaggerated so the space can be seen; in the mouth the margin lies close against the tooth. No fluid flow or composition is modeled.'],
 sulcular:['Sulcular epithelium','Non-keratinized epithelium lining the soft-tissue wall of the sulcus or pocket, continuous with the oral epithelium at the margin.','A thin authored layer facing the fluid space from the margin to the sulcus base. Thickness is exaggerated.'],
 junctional:['Junctional epithelium','The epithelial collar attached to the tooth at the base of the sulcus or pocket; widest coronally, tapering apically.','An authored wedge whose length is half the distance from the attachment preset to the crest preset — a display proportion, not a measurement.'],
 margin:['Gingival margin','The rounded free edge of the gingiva, encircling the tooth like a cuff.','The margin is swept around the root so it reads as a continuous collar. Height follows the margin preset.'],
 oral:['Oral gingival epithelium','Keratinized epithelium covering the outer gingiva; its connective-tissue interface forms rete ridges.','The wavy interface on the section is an artistic cue, not counted or measured ridges.'],
 connective:['Gingival connective tissue','Fibrous, vascular connective tissue beneath the epithelia; supracrestal fibres insert into cementum between the junctional epithelium and the bone crest.','The section pattern is an artistic fibre cue; fibre groups, directions, vessels and inflammation are not modeled.'],
 pdl:['Periodontal ligament','Fibrous connective tissue joining root cementum to the alveolar socket.','A uniform authored band, widened for visibility. Fibre orientation, stress and perfusion are not calculated.'],
 bone:['Alveolar crest and bone','Bone of the tooth socket; the crest is its coronal edge.','A rounded crest and outer plate. Crest height follows the preset and is not a simulated remodeling result.'],
 cuff:['Gingival margin · retained band','The marginal gingiva, kept while the deeper soft tissue is cut away.','In Biofilm surface the sulcular lining and connective tissue are removed to expose the tooth-facing plaque; the band marks where the margin lies.']
} as const;
export type CreviceDetail=keyof typeof DETAILS;
export const CREVICE_DETAILS:Record<CreviceDetail,{title:string;is:string;shown:string}>=Object.fromEntries(Object.entries(DETAILS).map(([k,[title,is,shown]])=>[k,{title,is,shown}])) as Record<CreviceDetail,{title:string;is:string;shown:string}>;
// Gentle artistic unevenness on the fluid-facing side of the plaque film only;
// not measured thickness, microcolony structure or species positions.
function plaqueRelief(f:ReturnType<typeof frame>){
 return (x:number,y:number,z:number)=>{
  if(x<f.ts(y)+.012)return 0;
  const n=.5+.25*Math.sin(y*5.1+z*1.3)*Math.sin(z*4.7-y*1.9)+.25*Math.sin(y*11.3+z*6.1)*Math.sin(z*9.7+y*2.3);
  return .045*n*ss(z,-HALF,-HALF+.35)*(1-ss(z,HALF-.35,HALF))*ss(y,f.plaqueBottom,f.plaqueBottom+.5)*(1-ss(y,f.plaqueTop-.4,f.plaqueTop));
 };
}
export function buildCrevice(health:Health):CrevicePart[]{
 const profiles=creviceProfiles(health),f=frame(health),all:CreviceView[]=['section','entrance','wall'],open:CreviceView[]=['section','entrance'];
 const tissue=(key:string,id:PocketId,name:string,color:string,finish:CreviceFinish,detail:CreviceDetail,views=all):CrevicePart=>({id,name,geometry:loft(profiles[key],key==='plaque'?plaqueRelief(f):undefined),color,finish,detail,views});
 const out:CrevicePart[]=[
  tissue('enamel','enamel','Enamel · coronal context','#e6dfcd','enamel','enamel'),
  tissue('dentin','dentin','Root dentin · tapering root','#e0c99c','hard','dentin'),
  tissue('pulp','pulp','Pulp chamber and canal · context','#b8655b','pulp','pulp'),
  tissue('cementum','cementum','Cementum · root surface','#cdb27b','hard','cementum'),
  tissue('pdl','pdl','Periodontal ligament · widened band','#bf7f74','ligament','pdl'),
  tissue('bone','bone','Alveolar crest and outer plate','#e2d5ba','bone','bone'),
  tissue('connective','connective','Gingival connective tissue','#cc887b','soft','connective',open),
  tissue('sulcular','epithelium','Sulcular epithelium · pocket lining','#e6a596','lining','sulcular',open),
  tissue('junctional','epithelium','Junctional epithelium · attachment','#f0c4b0','lining','junctional'),
  tissue('gingiva','gingiva','Gingival margin and oral epithelium','#df978b','soft','margin',open),
  tissue('cuff','gingiva','Gingival margin · retained band','#df978b','soft','cuff',['wall']),
  tissue('plaque','plaque','Attached plaque · illustrative film','#a99467','film','plaque')
 ];
 const s=cellScale(health);
 for(const p of microbialPlacements(health)){
  const t=TAXA.find(t=>t.id===p.id)!,g=microbialObject(p.id);
  g.rotateX((TAXA.indexOf(t)%3-1)*.45);g.scale(s,s,s);g.rotateY(-p.z/RADIUS);g.translate(...p.position);
  out.push({id:'plaque',name:t.name+' · illustrative 3D form',geometry:g,color:t.color,taxon:t.id,source:MORPHOLOGY[t.id].source,finish:'cell',views:all});
 }
 return out;
}
export interface CreviceLabel {key:string;id:PocketId;detail?:CreviceDetail;taxon?:string;text:string;short:string;point:[number,number,number];essential:boolean}
// Label anchors lie on the deliberate section face or on surfaces each camera
// faces. They are computed from the same profile curves as the meshes.
export function creviceLabels(view:CreviceView,health:Health,taxon='all'):CreviceLabel[]{
 const f=frame(health),{M,A,J,C,PT,B,ts,S,E,O,G,sM,pdl}=f,deep=health==='periodontitis';
 const label=(key:string,id:PocketId,detail:CreviceDetail|undefined,text:string,short:string,x:number,y:number,z:number,essential=false):CreviceLabel=>({key,id,detail,text,short,point:crevicePoint(x,y,z),essential});
 if(view==='section'){
  const yl=M-.45*(M-A),ys=M-.66*(M-A),yj=A-.3*(A-J),yc=(J+PT)/2,yp=C-.42*(C-B);
  return [
   label('enamel','enamel','enamel','Enamel','Enamel',(dentinWall(1.4)+enamelOuter(1.4))/2,1.4,HALF),
   label('margin','gingiva','margin','Gingival margin','Margin',sM+.3,M+.24,HALF,true),
   label('oral','gingiva','oral','Oral gingival epithelium','Oral epith.',(O(M-2.2)+G(M-2.2))/2,M-2.2,HALF),
   label('plaque','plaque','plaque','Attached plaque','Plaque',ts(M-.2*(M-A))+f.film(M-.2*(M-A))/2,M-.2*(M-A),HALF,true),
   label('lumen','lumen','lumen',deep?'Pocket · fluid space':'Sulcus · fluid space',deep?'Pocket':'Sulcus',ts(yl)+creviceGap(yl,health)*.55,yl,HALF,true),
   label('sulcular','epithelium','sulcular','Sulcular epithelium','Sulcular',S(ys)+.11,ys,HALF),
   label('junctional','epithelium','junctional','Junctional epithelium','Junctional',ts(yj)+f.junctional(yj)*.45,yj,HALF,true),
   label('connective','connective','connective','Connective tissue','Connective',(E(A)+G(A))/2,A-.2,HALF),
   label('cementum','cementum','cementum','Cementum','Cementum',dentinWall(yc)+cementum(yc)/2,yc,HALF),
   label('pdl','pdl','pdl','Periodontal ligament','Ligament',ts(yp)+pdl/2,yp,HALF,true),
   label('bone','bone','bone','Alveolar crest','Crest',ts(C)+pdl+.36,C-.3,HALF,true),
   label('dentin','dentin','dentin','Root dentin','Dentin',(dentinWall(M-3)+canal(M-3))/2,M-3,HALF),
   label('pulp','pulp','pulp','Pulp','Pulp',(canal(M-3)+CORE)/2+.15,M-3,HALF)
  ];
 }
 if(view==='entrance')return [
  label('margin','gingiva','margin','Gingival margin · cuff','Margin',sM+.32,M+.33,.4,true),
  label('lumen','lumen','lumen',deep?'Pocket entrance':'Sulcus entrance','Entrance',ts(M)+.13,M-.02,-.6,true),
  label('depth','lumen','lumen',deep?'Pocket depth at the cut face':'Sulcus depth at the cut face','Depth',ts(M-.4*(M-A))+creviceGap(M-.4*(M-A),health)*.5,M-.4*(M-A),HALF),
  label('plaque','plaque','plaque','Plaque at the margin','Plaque',ts(M+.22)+.04,M+.22,1.6),
  label('enamel','enamel','enamel','Enamel crown','Enamel',enamelOuter(2.2)-.02,2.2,0,true),
  label('oral','gingiva','oral','Oral gingiva','Gingiva',O(M-1.1),M-1.1,1.4)
 ];
 const out=[
  label('cuff','gingiva','cuff','Gingival margin (retained)','Margin',O(M-.35*f.cuffDepth),M-.35*f.cuffDepth,-.6,true),
  label('plaque','plaque','plaque','Attached plaque · tooth side','Plaque',ts(M-Math.max(1.25,.2*(M-A)))+f.film(M-Math.max(1.25,.2*(M-A))),M-Math.max(1.25,.2*(M-A)),2.5,true),
  label('junctional','epithelium','junctional','Junctional attachment','Junctional',ts(A-.35*(A-J))+f.junctional(A-.35*(A-J)),A-.35*(A-J),.4,true),
  label('cementum','cementum','cementum','Root surface · cementum','Cementum',ts((J+PT)/2),(J+PT)/2,-1.4),
  label('bone','bone','bone','Alveolar crest','Crest',ts(C)+pdl+.32,C,.8,true)
 ];
 if(taxon!=='all'){const p=microbialPlacements(health).find(p=>p.id===taxon),t=TAXA.find(t=>t.id===taxon);if(p&&t)out.push({key:'taxon',id:'plaque',taxon,text:t.short+' · illustrative form',short:t.short,point:p.position,essential:true});}
 return out;
}
export function creviceCamera(view:CreviceView,health:Health,taxon='all'){
 const f=frame(health),{M,C,B}=f,up=[0,1,0],healthy=health==='healthy';
 const pose=(target:number[],dir:number[],distance:number)=>{const l=Math.hypot(...dir);return {position:target.map((v,i)=>v+dir[i]/l*distance),target};};
 const basis=(z:number)=>{const r=creviceRadial(z);return {r,t:[-r[2],0,r[0]]};};
 const fit=(height:number)=>height/2/Math.tan(T.MathUtils.degToRad(17))*1.08;
 if(view==='entrance'){
  // Look down past the margin at the opened end: the sulcus shows as a V in the
  // cut face and the same narrow slit continues around the root as the cuff
  // curves away and fades.
  const z=HALF*.45,{r,t}=basis(z),target=crevicePoint(f.ts(M)+.4,M+1.1,z);
  return pose(target,r.map((v,i)=>v*.8+up[i]*1.55+t[i]*.95),healthy?18.5:22);
 }
 if(view==='wall'){
  const p=microbialPlacements(health).find(p=>p.id===taxon);
  if(p){const {r,t}=basis(p.z);return pose(p.position,r.map((v,i)=>v+up[i]*.04+t[i]*.1),healthy?7.5:15);}
  // The whole gallery between the retained margin and the crest.
  const top=M+1.1,bottom=C-.3,y=(top+bottom)/2,{r,t}=basis(0),target=crevicePoint(f.ts(y)+.2,y,0);
  // Level with the film so neither the margin band nor the junctional collar
  // hides the gallery rows in perspective.
  return pose(target,r.map((v,i)=>v+up[i]*.02+t[i]*.12),fit(top-bottom));
 }
 // Opened section: a three-quarter view of the deliberate section face, turned
 // so the outer gingiva and the swept cuff stay visible beside it.
 const top=5.3,bottom=B+.2,{r,t}=basis(HALF),target=crevicePoint(-1.75,(top+bottom)/2,HALF*.8);
 return pose(target,t.map((v,i)=>v+r[i]*.95+up[i]*.36),fit(top-bottom));
}
