import * as T from 'three';
import {mergeGeometries,mergeVertices} from 'three/addons/utils/BufferGeometryUtils.js';
import {TAXA,PRESETS,type Health} from './pocket-science';
import {MORPHOLOGY} from './pocket-morphology';
import type {PocketId} from './pocket-geometry';
export type CreviceView='section'|'entrance'|'wall';
export interface CrevicePart {id:PocketId;name:string;geometry:T.BufferGeometry;color:string;taxon?:string;source?:string;}
export const CREVICE_CLAIM='Magnified authored crevice; widths and cell sizes exaggerated independently. Example placements are not measured colonization, contacts or density-derived cells.';
export const toothWall=(y:number)=>-1.65+.34*Math.cos(y*.22)+.035*y;
// Shared curvature bends every interface through depth. These are authored
// display coordinates, not a reconstructed root or a physical-scale basis.
export function crevicePoint(x:number,y:number,z:number):[number,number,number]{
 const u=z/3.2,w=.93+.07*Math.cos(u*Math.PI/2);
 return [toothWall(y)+(x-toothWall(y))*w+.105*z*z+.045*z*Math.sin(y*.27),y+.22*Math.cos(z*.65)+.025*z*Math.sin(y*.24),z];
}
export function creviceLandmarks(health:Health){const p=PRESETS[health];return {margin:p.margin*1.6,attachment:p.attachment*1.6,crest:p.crest*1.6};}
export function creviceGap(y:number,health:Health){
 const {margin,attachment}=creviceLandmarks(health),t=T.MathUtils.clamp((margin-y)/(margin-attachment),0,1);
 return .16+(health==='healthy'?1.12:2.18)*Math.sin(t*Math.PI)**.85;
}
function loft(profile:T.Shape,depth=6.4){
 // Sample each curve by its parameter so shared tissue boundaries receive
 // exactly the same vertices even when their overall perimeters differ.
 const ring=profile.getPoints(36).slice(0,-1),steps=24;
 // Collinear samples along deliberate closure edges carry no silhouette
 // information and can make the cap triangulator emit zero-area triangles.
 for(let i=ring.length-1;i>=0;i--){const a=ring[(i+ring.length-1)%ring.length],b=ring[i],c=ring[(i+1)%ring.length];if(Math.abs((b.x-a.x)*(c.y-b.y)-(b.y-a.y)*(c.x-b.x))<1e-10)ring.splice(i,1);}
 const n=ring.length;
 if(T.ShapeUtils.isClockWise(ring))ring.reverse();
 const positions:number[]=[],uv:number[]=[],indices:number[]=[];
 for(let j=0;j<=steps;j++)for(let i=0;i<n;i++){
  const z=-depth/2+depth*j/steps,p=ring[i];positions.push(...crevicePoint(p.x,p.y,z));uv.push(i/n,j/steps);
 }
 for(let j=0;j<steps;j++)for(let i=0;i<n;i++){
  const a=j*n+i,b=j*n+(i+1)%n,c=(j+1)*n+i,d=(j+1)*n+(i+1)%n;
  indices.push(a,b,c,b,d,c);
 }
 const surfaceCount=indices.length,triangles=T.ShapeUtils.triangulateShape(ring,[]);
 // Separate cap vertices keep a clean normal break and a dedicated section
 // material. Their positions exactly match the loft's end rings.
 for(const end of [0,steps]){
  const offset=positions.length/3;
  for(const p of ring){positions.push(...crevicePoint(p.x,p.y,-depth/2+depth*end/steps));uv.push((p.x+6)/13,(p.y+16)/20);}
  for(const [a,b,c]of triangles)indices.push(...(end===0?[offset+a,offset+c,offset+b]:[offset+a,offset+b,offset+c]));
 }
 // Orient the entire closed volume consistently outward.
 let volume=0;for(let i=0;i<indices.length;i+=3){const a=indices[i]*3,b=indices[i+1]*3,c=indices[i+2]*3;volume+=positions[a]*(positions[b+1]*positions[c+2]-positions[b+2]*positions[c+1])+positions[a+1]*(positions[b+2]*positions[c]-positions[b]*positions[c+2])+positions[a+2]*(positions[b]*positions[c+1]-positions[b+1]*positions[c]);}
 if(volume<0)for(let i=0;i<indices.length;i+=3)[indices[i+1],indices[i+2]]=[indices[i+2],indices[i+1]];
 const g=new T.BufferGeometry();g.setAttribute('position',new T.Float32BufferAttribute(positions,3));g.setAttribute('uv',new T.Float32BufferAttribute(uv,2));g.setIndex(indices);g.addGroup(0,surfaceCount,0);g.addGroup(surfaceCount,indices.length-surfaceCount,1);g.computeVertexNormals();g.computeBoundingBox();g.userData.closedLoft=true;return g;
}
function ribbon(top:number,bottom:number,left:(y:number)=>number,right:(y:number)=>number,depth=6.4){
 const ys=Array.from({length:97},(_,i)=>top+(bottom-top)*i/96),s=new T.Shape();s.moveTo(left(top),top);for(const y of ys.slice(1))s.lineTo(left(y),y);for(const y of [...ys].reverse())s.lineTo(right(y),y);s.closePath();return loft(s,depth);
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
export function microbialPlacements(health:Health){
 const {margin,attachment}=creviceLandmarks(health),span=margin-attachment;
 // All identities use the same three display rows. No succession or depth niche
 // is assigned. Rows are a gallery arrangement, not positions from a specimen.
 return TAXA.map((t,i)=>{const row=Math.floor(i/3),y=margin-span*(.24+row*.27);return {id:t.id,position:crevicePoint(toothWall(y)+.62,y,-1.72+(i%3)*1.72)};});
}
export function buildCrevice(health:Health):CrevicePart[]{
 const {margin,attachment,crest}=creviceLandmarks(health),bottom=crest-3.3;
 const inner=(y:number)=>toothWall(y)+(y<attachment?.16:creviceGap(y,health));
 const lining=(y:number)=>inner(y)+.25;
 const coreInner=(y:number)=>y>=attachment?lining(y):toothWall(y)+.66;
 const outer=5.35,baseInner=coreInner(crest),floor=(x:number)=>crest-.24*Math.sin(Math.PI*(x-baseInner)/(outer-baseInner));
 const dentin=new T.Shape();dentin.moveTo(toothWall(3.2),3.2);
 for(let i=1;i<=120;i++){const y=3.2+(bottom-3.2)*i/120;dentin.lineTo(toothWall(y),y);}
 dentin.quadraticCurveTo(-3.6,bottom-.45,-3.85,bottom+.45);dentin.bezierCurveTo(-4.7,bottom+3,-5.35,-3,-4.55,2.75);dentin.quadraticCurveTo(-3.7,3.6,toothWall(3.2),3.2);dentin.closePath();
 const core=new T.Shape(),foldY=margin-.72;core.moveTo(lining(foldY),foldY);
 for(let i=1;i<=120;i++){const y=foldY+(crest-foldY)*i/120;core.lineTo(coreInner(y),y);}
 for(let i=1;i<=48;i++){const x=baseInner+(outer-baseInner)*i/48;core.lineTo(x,floor(x));}
 core.bezierCurveTo(6.1,crest+1.7,6.65,margin-3.15,4.0,margin-1.22);core.quadraticCurveTo(2.2,margin-.42,lining(foldY),foldY);core.closePath();
 // A continuous hood grows from the same core boundary. There are no sphere
 // ends or isolated tubes at the entrance. The top turns back into the core.
 const fold=new T.Shape();fold.moveTo(lining(foldY),foldY);
 for(let i=1;i<=32;i++){const y=foldY+(margin-foldY)*i/32;fold.lineTo(inner(y),y);}
 fold.quadraticCurveTo(inner(margin)+.12,margin+.29,inner(margin)+.68,margin+.13);
 fold.bezierCurveTo(1.45,margin-.05,3.1,margin-.42,4.0,margin-1.22);fold.quadraticCurveTo(2.2,margin-.42,lining(foldY),foldY);fold.closePath();
 const bone=new T.Shape();bone.moveTo(baseInner,crest);
 for(let i=1;i<=48;i++){const x=baseInner+(outer-baseInner)*i/48;bone.lineTo(x,floor(x));}
 bone.bezierCurveTo(5.6,crest-1.2,4.5,bottom+.3,2.4,bottom-.12);bone.quadraticCurveTo(.2,bottom-.4,toothWall(bottom)+.66,bottom);
 for(let i=1;i<=64;i++){const y=bottom+(crest-bottom)*i/64;bone.lineTo(toothWall(y)+.66,y);}bone.closePath();
 const out:CrevicePart[]=[
  {id:'dentin',name:'Root dentin · curved magnified section',geometry:loft(dentin),color:'#dcc5a0'},
  {id:'cementum',name:'Cementum · curved tooth-facing root surface',geometry:ribbon(0,bottom,toothWall,y=>toothWall(y)+.11),color:'#e9ddc2'},
  {id:'enamel',name:'Cervical enamel · coronal context',geometry:ribbon(3.2,0,toothWall,y=>toothWall(y)+.23),color:'#f2eee1'},
  {id:'connective',name:'Connective core · deliberate section faces',geometry:loft(core),color:'#c88376'},
  {id:'epithelium',name:'Sulcular lining and attachment · magnified',geometry:ribbon(foldY,attachment,inner,lining),color:'#e3a49a'},
  {id:'gingiva',name:'Continuous marginal gingival fold',geometry:loft(fold),color:'#c97870'},
  {id:'bone',name:'Alveolar crest and curved socket support',geometry:loft(bone),color:'#dbcfb8'},
  {id:'pdl',name:'Attachment-side ligament envelope · schematic',geometry:ribbon(attachment,bottom,y=>toothWall(y)+.16,y=>toothWall(y)+.66),color:'#cca399'}
 ];
 const plaqueYs=Array.from({length:65},(_,i)=>margin-.35+(attachment-margin+.55)*i/64);
 out.push({id:'plaque',name:'Attached plaque · illustrative coating',geometry:ribbon(plaqueYs[0],plaqueYs.at(-1)!,y=>toothWall(y)+.12,y=>toothWall(y)+.2,5.9),color:'#9b8e70'});
 for(const p of microbialPlacements(health)){
  const t=TAXA.find(t=>t.id===p.id)!,g=microbialObject(p.id);
  g.rotateX((TAXA.indexOf(t)%3-1)*.45);g.scale(health==='healthy'?.45:.85,health==='healthy'?.45:.85,health==='healthy'?.45:.85);g.translate(...p.position);
  out.push({id:'plaque',name:t.name+' · illustrative 3D form',geometry:g,color:t.color,taxon:t.id,source:MORPHOLOGY[t.id].source});
 }
 return out;
}
export function creviceCamera(view:CreviceView,health:Health,taxon='all'){
 const {margin,attachment,crest}=creviceLandmarks(health),middle=(margin+attachment)/2;
 if(view==='entrance')return {position:[10,margin+18,24],target:[.9,middle-.7,0]};
 if(view==='wall'){
  const p=microbialPlacements(health).find(p=>p.id===taxon)?.position??[-1.35,middle,0];
  return {position:[p[0]+(taxon==='all'?15:8),p[1]+1.25,p[2]+(taxon==='all'?2.3:2.8)],target:p};
 }
 const center=(margin+crest)/2;
 return {position:[12,center+3.0,health==='healthy'?25:36],target:[.65,center-.8,0]};
}
