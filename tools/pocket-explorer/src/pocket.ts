import * as T from 'three';
import {EffectComposer} from 'three/addons/postprocessing/EffectComposer.js';
import {RenderPass} from 'three/addons/postprocessing/RenderPass.js';
import {SSAOPass} from 'three/addons/postprocessing/SSAOPass.js';
import {OutputPass} from 'three/addons/postprocessing/OutputPass.js';
import {OrbitControls} from 'three/addons/controls/OrbitControls.js';
import {environment} from './materials';
import {pocketMaterials as materials,boneSectionTexture,connectiveSectionTexture} from './pocket-materials';
import {PocketViewCache,triangleCount,type PreparedPart} from './pocket-view-cache';
import {sceneLabelIds,visibleLabels,shortLabels,sceneCopy,additionalLabelTargets} from './pocket-layout';
import {buildOralSides,ANCHORS,PROFILES,type PocketId,type PocketPart} from './pocket-geometry';
import {TAXA,SOURCES,PRESETS,observation,validateSaved,gridRGBA,type Health,type Scale,type Layer,type Basis,type SavedInput} from './pocket-science';
import rawInput from './pocket-input.json';
import {REFERENCES,PARAMETERS,NERVE_STEPS,CROWN_REFERENCE} from './pocket-reference';
import {MORPHOLOGY,morphologyPanel} from './pocket-morphology';
import {buildCrevice,creviceCamera,creviceLandmarks,crevicePoint,toothWall,creviceGap,CREVICE_CLAIM,type CreviceView} from './pocket-crevice';
import {SEQUENCES,stageAt,sequenceSVG,sourceForStage,type Sequence} from './pocket-sequences';
import {neuralContext,FDI36_POSE,type ContextPart} from './pocket-context';
import {atlasContext,ATLAS_FDI36,ATLAS_SOURCE,ATLAS_CREDIT,type AtlasPart} from './pocket-atlas';

const get=<E extends HTMLElement=HTMLElement>(id:string)=>document.getElementById(id) as E;
const stage=get('stage'),canvas=get<HTMLCanvasElement>('scene'),labels=get('label-layer'),lens=get('lens');
stage.dataset.browser=navigator.userAgent;
const saved=validateSaved(rawInput as SavedInput);
const entryView=new URLSearchParams(location.search);
// The query override checks this same presentation policy without changing OS settings.
const reducedMotion=matchMedia('(prefers-reduced-motion: reduce)').matches||entryView.get('motion')==='reduce';
document.documentElement.dataset.reducedMotion=String(reducedMotion);
let health:Health=entryView.get('health')==='healthy'?'healthy':'periodontitis',scale:Scale=entryView.get('view')==='crevice'?'crevice':'pocket',layer:Layer='species',basis:Basis='evidence',taxon='all',selected:PocketId='plaque',cutaway=entryView.get('cutaway')!=='0',showLabels=true,inspector=true;
let renderer:T.WebGLRenderer|undefined,controls:OrbitControls|undefined,mats:ReturnType<typeof materials>|undefined,dirty=true,fallback=false;
const scene=new T.Scene(),camera=new T.PerspectiveCamera(34,1,.1,650),root=new T.Group();scene.add(root);
const connector=document.createElementNS('http://www.w3.org/2000/svg','svg');connector.classList.add('model-connector');connector.setAttribute('aria-hidden','true');stage.append(connector);
let boneCapMap:T.Texture|undefined,fiberCapMap:T.Texture|undefined;
let composer:EffectComposer|undefined,ao:SSAOPass|undefined;
const cache=new Map<string,PocketPart[]>(),display:T.Mesh[]=[],capMats=new Map<string,T.MeshPhysicalMaterial>();
const tissueAnchors=new Map<PocketId,T.Vector3>();
const remembered=new Map<Scale,{position:number[];target:number[]}>();
let creviceView:CreviceView='section';
let cameraTravel:{start:number;from:T.Vector3;to:T.Vector3;fromTarget:T.Vector3;toTarget:T.Vector3}|undefined;
function stopCameraTravel(){cameraTravel=undefined;stage.dataset.cameraTransition='idle';}
let neuro=false,sequence:Sequence='assembly',storyProgress=0,playing=false,storyLast=0;
const referenceMatrix=new T.Matrix4().makeScale(...CROWN_REFERENCE.scale);
const appearanceMats=new Map<string,T.MeshPhysicalMaterial>();
function appearanceMaterial(id:NonNullable<ContextPart['appearance']>){
 const facial=['skin','eye','iris'].includes(id);
 const key=facial?id+':'+(cutaway||neuro):id;
 if(!appearanceMats.has(key)){
  const color={skin:'#c39c81',lip:'#aa6760',eye:'#e1dfd0',iris:'#413d2e',vein:'#527eac',artery:'#bd4e44',nerve:'#d9ba70'}[id];
  const windowed=facial&&(cutaway||neuro);
  const m=new T.MeshPhysicalMaterial({color,roughness:id==='eye'?.24:id==='skin'?.61:.45,clearcoat:id==='lip'?.22:.06,side:T.DoubleSide,clippingPlanes:windowed?[new T.Plane(new T.Vector3(-1,0,0),0)]:[],transparent:false,opacity:1,depthWrite:true});
  if(id==='skin'||id==='lip'){m.map=id==='lip'?mats?.gingiva.map??null:null;m.bumpMap=mats?.gingiva.bumpMap??null;m.bumpScale=.012;m.color.set(color);}
  if(['nerve','artery','vein'].includes(id)){m.emissive.set(color);m.emissiveIntensity=.08;}
  appearanceMats.set(key,m);
 }return appearanceMats.get(key)!;
}
let measuring=false,measureStart=0,frameCount=0,renderMS:number[]=[];
const rank:Record<PocketId,number>={cheek:0,supragingival:8,bone:1,pdl:2,cementum:3,gingiva:4,connective:5,lumen:6,epithelium:7,plaque:8,enamel:9,dentin:10,pulp:11,vessels:12,nerve:13,palate:1,tongue:1};
const names:Record<PocketId,string>={supragingival:'Supragingival plaque',cheek:'Cheek side (buccal)',enamel:'Enamel',dentin:'Dentin',pulp:'Pulp & canals',cementum:'Cementum',pdl:'Periodontal ligament',gingiva:'Gingiva',bone:'Alveolar bone',vessels:'Vascular supply',nerve:'Neural supply',palate:'Palate',tongue:'Tongue side (lingual)',plaque:'Subgingival plaque',lumen:'Pocket fluid',epithelium:'Epithelial lining',connective:'Connective tissue'};
const descriptions:Record<PocketId,[string,string]>={
 enamel:['The mineralized outer crown; mature erupted enamel is acellular.','The crown shell and its cut edge are authored shapes, not a measured enamel-thickness map.'],
 dentin:['Mineralized tissue surrounding the pulp; odontoblast cell bodies lie on the pulp-facing boundary.','Coronal dentin and two root bodies are shown. Dentinal tubules are not resolved.'],
 pulp:['Living connective tissue containing vascular, neural and immune-associated populations.','Connected authored chamber, pulp horns and tapered canals. Canal number and geometry vary between real teeth.'],
 cementum:['The mineralized root covering at the attachment interface.','A thin shell around the roots; its thickness is illustrative.'],
 pdl:['Connective tissue between cementum and the alveolar socket.','A separate root envelope. Fiber directions, stress and perfusion are not calculated.'],
 gingiva:['Gum tissue surrounding the cervical tooth region.','A narrow entry follows the tooth-facing margin. The margin preset is independent from attachment and bone-crest presets.'],
 bone:['The alveolar socket supports the roots through the periodontal ligament.','Selected-site and opposing bone plus interradicular support. The crest is a preset, not a simulated remodeling result.'],
 vessels:['Blood vessels supply living pulp and periodontal tissues.','Schematic routes only; diameters, flow and oxygen delivery are not measured.'],
 nerve:['Neural routes provide context for tooth sensation.','Illustrative pathways; no neural activity is simulated.'],
 palate:['Hard/soft palate context forms the roof of the mouth.','Procedural oral-cavity context, not reconstructed tissue.'],
 tongue:['The tongue lies on the lingual (inner) side of the lower molars; its lateral border rests near the lingual cusps.','Translucent authored context in Tooth view: the lateral tongue border beside the lingual cusps. Shape and position are illustrative; tongue movement, saliva flow and papillae are not simulated.'],
 cheek:['For a lower first molar the outer (facial) side faces the cheek, called buccal. Lips border the front teeth, not this tooth.','Translucent authored context in Tooth view, drawn back over the mesial half: the cheek lining folds down into the buccal vestibule below the gum. Position and thickness are illustrative, with no measured registration.'],
 supragingival:['Plaque on the crown surface above the gingival margin, bathed by saliva.','A thin authored film over the cervical crown on the cheek and tongue sides and between them. It marks a compartment only: thickness and extent are illustrative, and no species positions are assigned to it.'],
 plaque:['Tooth-attached microbial biofilm below the gingival margin, distinct from the surrounding fluid lumen.','A thin film lines the sulcus on the cheek and tongue sides; the selected distal site keeps its own pocket ribbon. These are compartments, not species positions. Open biofilm detail for source-scoped observations across plaque thickness, or the separate assumed model grid.'],
 lumen:['The fluid space between the tooth-facing surface and the epithelial wall.','A slit-like opening and an enlarged illustrative interior. Widths are exaggerated so the space can be seen: in the mouth the gum lies close against the tooth and no open gap is visible. The true width is unknown here; it is not filled with stretched microbial layers.'],
 epithelium:['The tooth-associated epithelial boundary; the junctional attachment is distinct from a clinical probing endpoint.','Sulcular lining and an apical attachment region are shown as a separate tissue. Boundaries are not measured histology.'],
 connective:['Vascular connective tissue beneath the epithelial boundary.','A distinct cut-face region behind the lining. No immune migration, inflammation dynamics or tissue oxygen field is computed.']
};
function status(text:string){get('status').textContent=text;}
function remember(){if(controls&&scale!=='biofilm')remembered.set(scale,{position:camera.position.toArray(),target:controls.target.toArray()});}
function view(reset=false){
 if(!controls)return;
 stopCameraTravel();
 const stored=!reset&&remembered.get(scale);if(stored){camera.position.fromArray(stored.position);controls.target.fromArray(stored.target);}
 else if(scale==='face'){camera.position.set(170,74,400);controls.target.set(0,61,-10);}
 else if(scale==='mouth'){camera.position.set(95,24,265);controls.target.set(0,-3,0);}
 else if(scale==='crevice'){const pose=creviceCamera(creviceView,health,taxon);controls.target.fromArray(pose.target);camera.position.fromArray(pose.position).sub(controls.target).multiplyScalar(Math.max(1,1.12/camera.aspect)).add(controls.target);}
 else{
  const box=frameBox(),center=box.getCenter(new T.Vector3());
  const direction=new T.Vector3(.48,.25,1).normalize(),right=new T.Vector3().crossVectors(camera.up,direction).normalize(),up=new T.Vector3().crossVectors(direction,right).normalize();
  const tangent=Math.tan(T.MathUtils.degToRad(camera.fov/2)),fraction=stage.clientWidth<450?.43:.66;let distance=20;
  const point=new T.Vector3();
  for(const mesh of display){if(mesh.userData.context)continue;const p=mesh.geometry.getAttribute('position');for(let i=0;i<p.count;i++){
   point.fromBufferAttribute(p,i);if(!box.containsPoint(point))continue;point.sub(center);distance=Math.max(distance,Math.abs(point.dot(up))/(tangent*.92)+point.dot(direction),Math.abs(point.dot(right))/(tangent*camera.aspect*fraction)+point.dot(direction));
  }}
  controls.target.copy(center);camera.position.copy(center).addScaledVector(direction,distance*(scale==='tooth'?1.10:1));
 }
 controls.minDistance=scale==='crevice'?4:20;controls.update();dirty=true;updateLabels();
}
function travelCreviceCamera(){
 if(!controls)return;
 const from=camera.position.clone(),fromTarget=controls.target.clone();view(true);
 if(reducedMotion){stage.dataset.cameraTransition='reduced-motion';return;}
 const to=camera.position.clone(),toTarget=controls.target.clone();camera.position.copy(from);controls.target.copy(fromTarget);controls.update();
 cameraTravel={start:performance.now(),from,to,fromTarget,toTarget};stage.dataset.cameraTransition='moving';dirty=true;
}
// Translucent oral-side context never drives framing or label placement.
function modelBox(){const box=new T.Box3();for(const mesh of display)if(!mesh.userData.context)box.expandByObject(mesh);return box;}
// Pocket view frames the selected 36 distal pocket; neighbours run out of frame,
// which reads as continuous tissue rather than a mounted specimen.
const POCKET_FOCUS=new T.Box3(new T.Vector3(-0.6,-9.6,-5.4),new T.Vector3(8.4,5.6,5.8));
function frameBox(){if(scale!=='pocket')return modelBox();const box=modelBox();return box.isEmpty()?box:box.intersect(POCKET_FOCUS.clone().applyMatrix4(referenceMatrix));}
function material(id:PocketId){
 if(id in mats!)return mats![id as keyof typeof mats];
 const colors={plaque:'#b9a46a',supragingival:'#d4c79a',lumen:'#2c5862',epithelium:'#e9a397',connective:'#dba196',cheek:'#c98a8a'};
 if(id==='cheek')return new T.MeshPhysicalMaterial({color:colors.cheek,roughness:.55,clearcoat:.15,vertexColors:true,transparent:true,depthWrite:false,side:T.DoubleSide});
 return new T.MeshPhysicalMaterial({color:colors[id as keyof typeof colors],roughness:id==='lumen'?.3:.6,clearcoat:id==='lumen'?.3:.1,side:T.DoubleSide});
}
const extraMats=new Map<PocketId,T.MeshPhysicalMaterial>();
// Per-vertex concavity of a non-indexed surface, 0 (flat or convex) to 0.32 (groove).
function concavity(g:T.BufferGeometry){
 const p=g.getAttribute('position'),n=g.getAttribute('normal'),id=new Map<string,number>(),of:number[]=[];
 for(let i=0;i<p.count;i++){const k=`${Math.round(p.getX(i)*1e4)},${Math.round(p.getY(i)*1e4)},${Math.round(p.getZ(i)*1e4)}`;if(!id.has(k))id.set(k,id.size);of.push(id.get(k)!);}
 const sum=new Float32Array(id.size*3),nor=new Float32Array(id.size*3),pos=new Float32Array(id.size*3),cnt=new Float32Array(id.size),len=new Float32Array(id.size);
 for(let i=0;i<p.count;i++){const v=of[i];pos.set([p.getX(i),p.getY(i),p.getZ(i)],v*3);nor[v*3]+=n.getX(i);nor[v*3+1]+=n.getY(i);nor[v*3+2]+=n.getZ(i);}
 for(let t=0;t<p.count;t+=3)for(let a=0;a<3;a++){const v=of[t+a];for(const b of [(a+1)%3,(a+2)%3]){const w=of[t+b];for(let c=0;c<3;c++)sum[v*3+c]+=pos[w*3+c];cnt[v]++;len[v]+=Math.hypot(pos[w*3]-pos[v*3],pos[w*3+1]-pos[v*3+1],pos[w*3+2]-pos[v*3+2]);}}
 const out=new Float32Array(p.count);
 for(let i=0;i<p.count;i++){const v=of[i];if(!cnt[v])continue;const nl=Math.hypot(nor[v*3],nor[v*3+1],nor[v*3+2])||1;let d=0;for(let c=0;c<3;c++)d+=(sum[v*3+c]/cnt[v]-pos[v*3+c])*nor[v*3+c]/nl;out[i]=Math.min(.32,Math.max(0,d/(len[v]/cnt[v])*1.6));}
 return out;
}
let tongueContext:T.MeshPhysicalMaterial|undefined;
function tongueMat(){return tongueContext??=new T.MeshPhysicalMaterial({color:'#c4777b',roughness:.5,clearcoat:.2,vertexColors:true,transparent:true,depthWrite:false,side:T.DoubleSide});}
function mat(id:PocketId){if(id in mats!)return material(id);if(!extraMats.has(id))extraMats.set(id,material(id));return extraMats.get(id)!;}
function clear(){for(const mesh of display){root.remove(mesh);if(mesh.userData.ephemeral){mesh.geometry.dispose();if(!mesh.userData.appearance)for(const m of Array.isArray(mesh.material)?mesh.material:[mesh.material])m.dispose();}}display.length=0;}
function anchorTarget(id:PocketId):[number,number,number]{
 if(id==='bone'&&!cutaway)return new T.Vector3(-6.4,-8.0,4.6).applyMatrix4(referenceMatrix).toArray() as [number,number,number];
 const anchor=[...(ANCHORS[id]??additionalLabelTargets[id]!)] as [number,number,number];if(['plaque','lumen','epithelium'].includes(id)){
   anchor[1]=id==='epithelium'?PRESETS[health].attachment+.10:(PRESETS[health].margin+PRESETS[health].attachment)/2;
   const y=anchor[1],t=Math.max(0,Math.min(1,(-y-2)/11.60));
   const neck=y>=-4.4?4.65*Math.sqrt(Math.max(0,1-((y+1.2)/3.2)**2)):0;
   const rootX=y<=-2?2+1.05*Math.sin(t*Math.PI/2)+.85*t*t+1.74*(1-.5*t)*(1-t**3.2)**.42+.03:0;
   anchor[0]=Math.max(neck,rootX)+.13+(id==='plaque'?.035:id==='lumen'?.24:.44);anchor[2]=0;
 }
 return new T.Vector3(...anchor).applyMatrix4(referenceMatrix).toArray() as [number,number,number];
}
function decorate(part:PreparedPart){
 const g=part.geometry,cap=part.cap;
 // Presentation only: original canonical meshes and scientific fields stay intact.
 g.applyMatrix4(referenceMatrix);cap?.applyMatrix4(referenceMatrix);
  if(part.id==='enamel'){
   // Optical cervical warmth, authored art rather than a measured mineral map.
   // Fissures and fossae are darkened slightly so grooves read at a glance.
   for(const mesh of [g,cap].filter(Boolean) as T.BufferGeometry[]){const p=mesh.getAttribute('position'),cavity=mesh===g?concavity(mesh):null,colors=[];for(let i=0;i<p.count;i++){const tint=new T.Color('#ddd0b4').lerp(new T.Color('#f2eee6'),T.MathUtils.clamp((p.getY(i)+.3)/5.2,0,1)**.8);if(cavity)tint.multiplyScalar(1-cavity[i]);colors.push(tint.r,tint.g,tint.b);}mesh.setAttribute('color',new T.Float32BufferAttribute(colors,3));}
   mat('enamel').vertexColors=true;
  }
 if(part.id==='gingiva'){
  const p=g.getAttribute('position'),colors:number[]=[],attached:number[]=[];
  for(let i=0;i<p.count;i++){
   const y=p.getY(i),plate=Math.abs(p.getZ(i)-.006*p.getX(i)**2)>4.35;
   const t=plate?T.MathUtils.smoothstep(y,-5.4,-4.25):1;
   const tint=new T.Color('#e6b2ae').lerp(new T.Color('#fff1e8'),t);
   colors.push(tint.r,tint.g,tint.b);attached.push(t);
  }
  g.setAttribute('color',new T.Float32BufferAttribute(colors,3));
  g.setAttribute('attachedGingiva',new T.Float32BufferAttribute(attached,1));
 }
}
const pocketCache=new PocketViewCache(decorate);
const anchorCache=new Map<string,Map<PocketId,T.Vector3>>();
function geometry(){
 updateFallback();if(!renderer||!mats)return;
 const started=performance.now();clear();
 if(scale==='crevice'){
  for(const part of buildCrevice(health)){
   const muted=part.taxon&&taxon!=='all'&&part.taxon!==taxon;
   const soft=['epithelium','connective','gingiva','pdl'].includes(part.id);
   const m=new T.MeshPhysicalMaterial({color:muted?'#7e8c87':part.color,roughness:part.taxon?.46:soft?.48:part.id==='enamel'?.3:.57,clearcoat:part.taxon?.12:soft?.16:.06,clearcoatRoughness:.55,sheen:soft?.22:0,sheenColor:'#d88d80',sheenRoughness:.85,specularIntensity:.65,envMapIntensity:.22,side:T.FrontSide});
   if(part.taxon){m.emissive.set(muted?'#000000':part.color);m.emissiveIntensity=.025;}
   if(soft){m.bumpMap=mats.gingiva.bumpMap;m.bumpScale=.014;}
   const section=m.clone();section.roughness=.84;section.clearcoat=0;section.sheen=0;section.specularIntensity=.25;section.envMapIntensity=.08;section.bumpScale=.008;
   // Deliberate section faces use restrained artistic tissue cues. These
   // patterns do not resolve measured fibers, trabeculae or cell-scale anatomy.
   if(['connective','bone'].includes(part.id)){
    section.map=part.id==='bone'?boneCapMap!:fiberCapMap!;section.color.set(part.id==='bone'?'#e9ddc2':'#de9b8b');section.bumpMap=null;
    const uv=part.geometry.getAttribute('uv');for(let i=0;i<uv.count;i++)uv.setXY(i,uv.getX(i)*4,uv.getY(i)*4);
    section.customProgramCacheKey=()=> 'marse-crevice-section-cues-1';
    section.onBeforeCompile=shader=>{shader.fragmentShader=shader.fragmentShader.replace('#include <map_fragment>',T.ShaderChunk.map_fragment.replace('diffuseColor *= sampledDiffuseColor;','diffuseColor *= mix(vec4(1.0), sampledDiffuseColor, 0.28);'));};
   }
   const mesh=new T.Mesh(part.geometry,part.taxon?m:[m,section]);if(part.taxon)section.dispose();mesh.name=part.name;mesh.userData={tissue:part.id,taxon:part.taxon,morphologySource:part.source,ephemeral:true,crevice:true};mesh.castShadow=true;mesh.receiveShadow=true;root.add(mesh);display.push(mesh);
  }
  stage.dataset.creviceSelection=health+':'+taxon;stage.dataset.microbialObjects=JSON.stringify(display.filter(m=>m.userData.taxon).map(m=>({taxon:m.userData.taxon,name:m.name,source:m.userData.morphologySource})));stage.dataset.creviceClaim=CREVICE_CLAIM;stage.dataset.triangles=String(display.reduce((n,m)=>n+triangleCount(m.geometry),0));stage.dataset.openContours='0';stage.dataset.closedContours='9 authored closed tissue lofts; software topology checks only';stage.dataset.geometryMs=(performance.now()-started).toFixed(2);stage.dataset.geometryCacheHit='false';dirty=true;updateLabels();return;
 }
 let open=0,loops=0,triangles=0;
 let parts:PreparedPart[],hit=false;
 if(scale==='mouth'||scale==='face'){
  const key=`oral:${health}:${scale}:${cutaway}`;hit=cache.has(key);
  if(!hit){
   cache.set(key,atlasContext(scale));
  }
  parts=cache.get(key)!;
 }else{
  hit=pocketCache.has(health,cutaway);const ready=pocketCache.get(health,cutaway);
  parts=ready.parts;open=ready.openChains;loops=ready.loops;
 }
 if(scale==='tooth'&&!cache.has('sides'))cache.set('sides',buildOralSides().map(p=>({...p,geometry:p.geometry.applyMatrix4(referenceMatrix)})));
 const context:ContextPart[]=scale==='tooth'?cache.get('sides')!:[];
 const routes=neuro&&scale!=='mouth'&&scale!=='face'?neuralContext(false):[];
 if(neuro&&scale!=='mouth'&&scale!=='face')for(const p of routes)p.geometry.applyMatrix4(referenceMatrix);
 for(const part of [...parts,...context,...routes] as (PreparedPart&ContextPart&Partial<AtlasPart>)[]){const isContext=context.includes(part);
  const g=part.geometry,cap=part.cap;
  // Shared enamel and gingiva materials require attributes on context surfaces too.
  if(!g.getAttribute('color'))g.setAttribute('color',new T.Float32BufferAttribute(new Float32Array(g.getAttribute('position').count*3).fill(1),3));
  if(part.id==='gingiva'&&!g.getAttribute('attachedGingiva'))g.setAttribute('attachedGingiva',new T.Float32BufferAttribute(new Float32Array(g.getAttribute('position').count).fill(1),1));
  if(neuro&&part.name.startsWith('Schematic ')&&['nerve','vessels'].includes(part.id))continue;
  const appearance=part.appearance;
  const styled=appearance?appearanceMaterial(appearance):isContext&&part.id==='tongue'?tongueMat():mat(part.id);
  const mesh=new T.Mesh(g,styled);mesh.name=part.name;mesh.userData={tissue:part.id,context:isContext,ephemeral:routes.includes(part),appearance,atlasIdentity:part.atlasIdentity};mesh.castShadow=!isContext;mesh.receiveShadow=true;if(isContext)mesh.renderOrder=30;root.add(mesh);display.push(mesh);triangles+=triangleCount(g);
  if(cap){if(cap.getAttribute('position').count){
   const capKey=part.id==='bone'?part.name:part.id;
   if(!capMats.has(capKey)){
    const m=mat(part.id).clone();
    m.map=part.id==='bone'?(part.name.includes('Trabecular')?boneCapMap!:null):part.id==='connective'?fiberCapMap!:mat(part.id).map;
    if(['gingiva','connective'].includes(part.id))m.color.set('#ffffff');
    // Cut faces are matte and slightly lighter than their glossy exteriors.
    m.bumpMap=part.id==='bone'?m.map:null;m.bumpScale=.04;m.clearcoat=0;m.sheen=0;m.roughness=.86;m.specularIntensity=.3;
    if(part.id==='dentin')m.color.set('#f3e2c6');if(part.id==='bone')m.color.set(part.name.includes('Trabecular')?'#ffffff':'#d8cfbd');if(part.id==='pulp')m.color.set('#f0d4cc');if(part.id==='enamel')m.color.set('#f4eee2');m.depthWrite=true;m.envMapIntensity=.06;m.polygonOffset=true;m.polygonOffsetFactor=-rank[part.id];m.polygonOffsetUnits=-rank[part.id];m.userData.tissue=part.id;capMats.set(capKey,m);
   }
   const surface=new T.Mesh(cap,capMats.get(capKey));surface.name=part.name+' section';surface.userData={tissue:part.id};surface.renderOrder=10+rank[part.id];root.add(surface);display.push(surface);triangles+=cap.getAttribute('position').count/3;
  }}
 }
 if(scale==='mouth'||scale==='face'){
  const ring=new T.Mesh(new T.TorusGeometry(3.8,.18,12,64),new T.MeshBasicMaterial({color:'#c8dcba'}));ring.rotation.x=-Math.PI/2;ring.position.copy(ATLAS_FDI36);ring.userData={tissue:'plaque',focus:true,ephemeral:true};root.add(ring);display.push(ring);
 }
 tissueAnchors.clear();const anchorKey=(scale==='mouth'||scale==='face')?`oral:${health}:${scale}:${neuro}`:`${health}:${cutaway}:${scale}`;
 const knownAnchors=anchorCache.get(anchorKey);
 if(knownAnchors){for(const [id,p]of knownAnchors)tissueAnchors.set(id,p);}
 else {
 const anchorCosts=new Map<PocketId,number>();
 for(const mesh of display){const id=mesh.userData.tissue as PocketId;if(!ANCHORS[id]&&!additionalLabelTargets[id])continue;const target=new T.Vector3(...anchorTarget(id)),p=mesh.geometry.getAttribute('position'),point=new T.Vector3();
  for(let i=0;i<p.count;i++){point.fromBufferAttribute(p,i);const d=point.distanceToSquared(target);if(d<(anchorCosts.get(id)??Infinity)){anchorCosts.set(id,d);tissueAnchors.set(id,point.clone());}}
 }
 anchorCache.set(anchorKey,new Map(tissueAnchors));
 }
 stage.dataset.geometryMs=(performance.now()-started).toFixed(2);
 stage.dataset.geometryCacheHit=String(hit);stage.dataset.geometryBuilds=String(pocketCache.builds);
 stage.dataset.tissueAnchors=JSON.stringify(Object.fromEntries([...tissueAnchors].map(([id,p])=>[id,p.toArray()])));
 stage.dataset.openContours=scale==='mouth'||scale==='face'?'source mesh; not sectioned':String(open);stage.dataset.closedContours=scale==='mouth'||scale==='face'?'not sectioned':String(loops);stage.dataset.triangles=String(Math.round(triangles));
 highlight();dirty=true;updateLabels();
 if(open)status('Some section contours are open; geometry requires correction.');
}
function highlight(){
 if(!mats)return;
 if(scale==='crevice'){for(const mesh of display){mesh.visible=creviceView!=='wall'||!['connective','epithelium'].includes(mesh.userData.tissue);if(!mesh.userData.taxon)for(const material of Array.isArray(mesh.material)?mesh.material:[mesh.material]){const m=material as T.MeshPhysicalMaterial;m.emissive.set(mesh.userData.tissue===selected&&layer==='structure'?'#71442a':'#000000');m.emissiveIntensity=.1;}}dirty=true;return;}
 // One tissue at a time: the selected tissue and its cut faces glow warmly.
 for(const id of Object.keys(rank) as PocketId[]){const m=mat(id);m.emissive.set(id===selected?'#5a3612':'#000000');m.emissiveIntensity=.32;for(const cap of capMats.values())if(cap.userData.tissue===id){cap.emissive.copy(m.emissive);cap.emissiveIntensity=.30;}}
 if(layer==='species'){mat('plaque').emissive.set('#a17c36');mat('plaque').emissiveIntensity=.28;}
 dirty=true;
}
function updateFallback(){
 const v=PRESETS[health];get('fallback-lumen').setAttribute('d',PROFILES[health].lumen);get('fallback-gum').setAttribute('d',PROFILES[health].gingiva);
 get('fallback-bone').setAttribute('d',`M10.2,${v.crest} L5.05,${v.crest} Q5.6,-11 5,-14.5 L2,-16.2 L10.3,-16.2 Z`);
 get('fallback-plaque').setAttribute('d',`M4.45,${v.margin-.2} L4.55,${v.margin-.2} L4.6,${v.attachment} L4.48,${v.attachment} Z`);
}
const labelIds=sceneLabelIds;
const labelElements=new Map<PocketId,HTMLButtonElement>();
const leaders=document.createElementNS('http://www.w3.org/2000/svg','svg');leaders.setAttribute('class','leaders');leaders.setAttribute('aria-hidden','true');labels.append(leaders);
for(const id of labelIds){const b=document.createElement('button');b.textContent=names[id];b.setAttribute('aria-label','Inspect '+names[id]);b.addEventListener('click',()=>selectTissue(id));labels.append(b);labelElements.set(id,b);}
function updateLabels(){
 camera.updateMatrixWorld();
 if(controls){stage.dataset.camera=JSON.stringify(camera.position.toArray());stage.dataset.cameraTarget=JSON.stringify(controls.target.toArray());}
 connector.replaceChildren();const connected=basis==='model'&&layer==='species'&&scale!=='mouth'&&scale!=='face'&&scale!=='biofilm'&&scale!=='crevice';connector.toggleAttribute('hidden',!connected);
 if(connected){
  const p=new T.Vector3(4.55,(PRESETS[health].margin+PRESETS[health].attachment)/2,0).applyMatrix4(referenceMatrix).project(camera),r=lens.getBoundingClientRect(),s=stage.getBoundingClientRect(),x=fallback?stage.clientWidth*.64:(p.x+1)/2*stage.clientWidth,y=fallback?stage.clientHeight*.45:(-p.y+1)/2*canvas.clientHeight+62;
  if(Number.isFinite(x)&&Number.isFinite(y)){  connector.setAttribute('viewBox',`0 0 ${stage.clientWidth} ${stage.clientHeight}`);const line=document.createElementNS('http://www.w3.org/2000/svg','path');line.setAttribute('d',`M${x},${y} L${r.left-s.left+9},${r.top-s.top+9}`);line.setAttribute('stroke','#a4bda8');line.setAttribute('stroke-opacity','.4');line.setAttribute('stroke-dasharray','3 5');line.setAttribute('fill','none');connector.append(line);}
 }
 leaders.replaceChildren();
 labels.hidden=!showLabels||scale==='biofilm';
 if(scale==='crevice'){
  const {margin,attachment,crest}=creviceLandmarks(health),middle=(margin+attachment)/2;
  const points:Partial<Record<PocketId,number[]>>={cementum:crevicePoint(toothWall(middle)+.08,middle,3.2),lumen:crevicePoint(toothWall(middle)+creviceGap(middle,health)/2,middle-1.1,3.2),epithelium:crevicePoint(toothWall(middle)+creviceGap(middle,health)+.14,middle,3.2),gingiva:crevicePoint(toothWall(margin)+.4,margin+.1,3.2),bone:crevicePoint(2.2,crest-.35,3.2),pdl:crevicePoint(toothWall((attachment+crest)/2)+.4,(attachment+crest)/2,3.2)};
  for(const [id,b]of labelElements){b.hidden=!points[id]||creviceView==='wall';b.setAttribute('aria-label','Inspect '+names[id]);if(b.hidden)continue;const p=new T.Vector3(...points[id] as [number,number,number]).project(camera);b.textContent=id==='cementum'?'Tooth surface':id==='epithelium'?'Tissue lining':id==='lumen'?'Fluid space':id==='gingiva'?'Pocket entrance':id==='pdl'?'Ligament':'Bone crest';const ax=(p.x+1)/2*stage.clientWidth,ay=(-p.y+1)/2*canvas.clientHeight+62,dx=id==='cementum'||id==='pdl'?-72:id==='epithelium'?88:id==='gingiva'?48:0,dy=id==='gingiva'?-28:id==='lumen'?35:0,x=T.MathUtils.clamp(ax+dx,58,stage.clientWidth-68),y=T.MathUtils.clamp(ay+dy,92,stage.clientHeight-136);b.style.left=x+'px';b.style.top=y+'px';leaders.setAttribute('viewBox',`0 0 ${stage.clientWidth} ${stage.clientHeight}`);const line=document.createElementNS('http://www.w3.org/2000/svg','path');line.setAttribute('d',`M${ax},${ay}L${x},${y}`);line.setAttribute('stroke','#b5c4b1');line.setAttribute('stroke-opacity','.45');line.setAttribute('fill','none');leaders.append(line);}
  return;
 }
 if(scale==='mouth'||scale==='face'){
  for(const [id,b]of labelElements){b.setAttribute('aria-label','Inspect '+names[id]);b.hidden=id!=='plaque';if(id==='plaque'){b.textContent='Selected FDI 36 →';b.setAttribute('aria-label','Open FDI 36 tooth detail');const p=ATLAS_FDI36.clone().project(camera);b.style.left=((p.x+1)/2*stage.clientWidth)+'px';b.style.top=((-p.y+1)/2*canvas.clientHeight+62+24)+'px';}}
  return;
 }
 const h=canvas.clientHeight||stage.clientHeight;
 const bounds=frameBox(),projected:Array<T.Vector3>=[];
 for(const x of [bounds.min.x,bounds.max.x])for(const y of [bounds.min.y,bounds.max.y])for(const z of [bounds.min.z,bounds.max.z])projected.push(new T.Vector3(x,y,z).project(camera));
 const xs=projected.map(p=>(p.x+1)/2*stage.clientWidth).filter(Number.isFinite),minX=xs.length?Math.min(...xs):0,maxX=xs.length?Math.max(...xs):stage.clientWidth;
 stage.dataset.modelBounds=JSON.stringify({left:minX,right:maxX});
 leaders.setAttribute('viewBox',`0 0 ${stage.clientWidth} ${stage.clientHeight}`);
 const compact=stage.clientWidth<560,visible=visibleLabels(stage.clientWidth,cutaway,scale,selected);
 stage.dataset.labelMode=compact?'compact':'full';get('labels').textContent=compact?'Labels · compact':'Labels';
 const left:Array<{b:HTMLButtonElement;x:number;y:number;ax:number;ay:number}>=[],right:Array<{b:HTMLButtonElement;x:number;y:number;ax:number;ay:number}>=[];
 for(const [id,b]of labelElements){b.setAttribute('aria-label','Inspect '+names[id]);b.hidden=!visible.includes(id);if(b.hidden)continue;b.textContent=(compact?shortLabels[id]:undefined)||({pdl:'Ligament',lumen:'Pocket fluid · width exaggerated',epithelium:'Junctional lining',supragingival:'Supragingival plaque',cheek:'Cheek side · buccal',tongue:'Tongue side · lingual'} as Partial<Record<PocketId,string>>)[id]||names[id];b.dataset.selected=String(id===selected);
  if(fallback){const positions:Record<string,[number,number]>={enamel:[23,22],pulp:[38,39],gingiva:[80,36],plaque:[61,45],epithelium:[79,55],pdl:[27,71],bone:[82,76],dentin:[20,32],cementum:[24,60],lumen:[80,47],supragingival:[33,27]};if(!positions[id]){b.hidden=true;continue;}const fh=stage.classList.contains('with-lens')?stage.clientHeight-247:stage.clientHeight-62;b.style.left=positions[id][0]+'%';b.style.top=(positions[id][1]*fh/100+62)+'px';continue;}
  const anchor=(tissueAnchors.get(id)||new T.Vector3(...anchorTarget(id))).toArray() as [number,number,number];
  const p=new T.Vector3(...anchor).project(camera),isLeft=['enamel','dentin','pulp','gingiva','supragingival','cementum','tongue'].includes(id)||(id==='bone'&&!cutaway);
  const ax=(p.x+1)/2*stage.clientWidth,ay=(-p.y+1)/2*h+62;
  if(!Number.isFinite(ax)||!Number.isFinite(ay)){b.hidden=true;continue;}
  const half=b.offsetWidth/2;const item={b,ax,ay,x:T.MathUtils.clamp(isLeft?minX-half-10:maxX+half+10,half+8,stage.clientWidth-half-8),y:T.MathUtils.clamp(ay,82,h+38)};(isLeft?left:right).push(item);
 }
 for(const list of [left,right]){list.sort((a,b)=>a.y-b.y);let last=60;for(let i=0;i<list.length;i++)list[i].y=Math.min(list[i].y,h+38-(list.length-1-i)*32);for(const item of list){item.y=Math.max(item.y,last+32);last=item.y;item.b.style.left=item.x+'px';item.b.style.top=item.y+'px';
  const line=document.createElementNS('http://www.w3.org/2000/svg','line');line.setAttribute('x1',String(item.ax));line.setAttribute('y1',String(item.ay));line.setAttribute('x2',String(item.x));line.setAttribute('y2',String(item.y));line.setAttribute('stroke','#9db9ac');line.setAttribute('stroke-opacity','.55');line.setAttribute('stroke-width','1');leaders.append(line);
  const dot=document.createElementNS('http://www.w3.org/2000/svg','circle');dot.setAttribute('cx',String(item.ax));dot.setAttribute('cy',String(item.ay));dot.setAttribute('r','2');dot.setAttribute('fill','#d8bc8b');leaders.append(dot);
 }}
}
function addSource(key:keyof typeof SOURCES){const s=SOURCES[key],a=document.createElement('a');a.href=s.url;a.textContent=s.title+' ↗';a.target='_blank';a.rel='noopener noreferrer';get('inspect-sources').append(a);}
function showInspector(){inspector=true;get('inspector').hidden=false;get('workspace').classList.remove('inspector-closed');get('show-inspector').hidden=true;resize();}
function inspect(title:string,is:string,shown:string,support:string,tag:string,source:keyof typeof SOURCES){
 get('inspect-title').textContent=title;get('inspect-is').textContent=is;get('inspect-shown').textContent=shown;get('inspect-support').textContent=support;get('inspect-tag').textContent=tag;get('inspect-sources').replaceChildren();addSource(source);
}
function inspectSpecies(){
 if(scale==='mouth'||scale==='face'){const part=display.find(m=>m.userData.atlasIdentity?.fdi===36);if(part)inspectAtlas(part.userData.atlasIdentity);return;}
 if(taxon==='all'){
  inspect('Attached biofilm','A structured microbial community attached to the tooth surface, distinct from the fluid space.',basis==='model'?'Nine taxon-density fields from a preserved teaching endpoint. Dominant colors do not mean the other taxa are absent.':'Qualitative organization across attached plaque. The schematic glyphs represent identity and region, not counts, sizes or contacts.',health==='healthy'?'The selected spatial study does not establish healthy-site positions. All species locations remain unresolved in this state.':SOURCES.plaque.method+' '+SOURCES.plaque.limit,basis==='model'?'Saved teaching fields · assumed grid':'Source-scoped observations · schematic detail','plaque');
 }else{
  const t=TAXA.find(t=>t.id===taxon)!,o=observation(taxon,health);
  inspect(t.name,`${t.level==='genus/group'?'Genus/group-level':'Species-level'} identity in the preserved teaching panel. This selected panel is not the complete oral microbiome.`,basis==='model'?'The original '+health+' density array in its unchanged 32 × 72 grid. Contrast is normalized to this taxon; colors cannot compare absolute abundance across selections.':o.text,basis==='model'?'The saved run uses imposed niche and inoculum assumptions. It supplies no measured cell positions or anatomical transformation.':SOURCES.plaque.method+' '+SOURCES.plaque.limit,basis==='model'?'Assumed model coordinates':o.resolved?'Human periodontitis observation':'Location unresolved','plaque');
 }
 if(basis==='model'){
  get('inspect-support').textContent='Preserved teaching output with imposed niche and inoculum assumptions. The spatial study supports the separate evidence view; it does not validate these densities or their anatomical positions.';
  get('inspect-sources').textContent=`Local source: ${saved.source.file} · SHA-256 ${saved.source.sha256}. ${saved.source.precision}.`;
 }
}
function inspectAtlas(identity:AtlasPart['atlasIdentity'],reveal=false){
 if(reveal)showInspector();inspect(identity.fdi?'FDI '+identity.fdi+' · atlas tooth':identity.name,'A named part of the BodyParts3D adult male reference atlas.','Source geometry in a shared coordinate system. Lower-jaw opening is an authored display pose; whole-tooth meshes do not separate enamel, dentin or pulp. The pocket detail is an independent specimen.','Source element '+identity.element+' · '+identity.concept+'. '+ATLAS_CREDIT+'. Reduced atlas meshes are not patient measurements or a population average. Expert review pending.','Reference atlas · authored opening pose','anatomy');get('inspect-sources').replaceChildren();const a=document.createElement('a');a.href=ATLAS_SOURCE;a.textContent='BodyParts3D source and identities ↗';a.target='_blank';a.rel='noopener noreferrer';get('inspect-sources').append(a);get('nerve-inspector').hidden=true;
}
function selectTissue(id:PocketId){
 if(scale==='crevice'){selected=id;layer='structure';showInspector();refresh();return;}
 if((scale==='mouth'||scale==='face')&&id==='plaque'){setScale('tooth');return;}
 selected=id;get<HTMLSelectElement>('tissue').value=id;showInspector();if(scale==='mouth'||scale==='face'){const m=display.find(m=>m.userData.atlasIdentity&&m.userData.tissue===id);if(m)inspectAtlas(m.userData.atlasIdentity);return;}if(id==='plaque'&&layer==='species'){inspectSpecies();if(scale!=='biofilm')setScale('biofilm');}
 else if(id==='supragingival')inspect(names[id],...descriptions[id],'This viewer assigns no species positions to supragingival plaque. The selected spatial observations concern attached subgingival biofilm in advanced periodontitis and are not transferred here.','Authored compartment · species locations unresolved','plaque');
 else inspect(names[id],...descriptions[id],SOURCES.anatomy.method+' '+SOURCES.anatomy.limit,'Authored anatomy · expert review pending','anatomy');
 get('nerve-inspector').hidden=!(id==='nerve'||id==='vessels');
 if(id==='nerve'||id==='vessels'){neuro=true;get('neuro').setAttribute('aria-pressed','true');geometry();renderNeuralInspector();}
 if(id==='enamel'){get('inspect-shown').textContent='The five-cusp crown is reference-sized to 11.13 mm mesiodistally and 10.47 mm buccolingually in its local frame. Remaining tissue dimensions are authored. This is a representative composite, not a measured tooth.';get('inspect-sources').replaceChildren();const a=document.createElement('a');a.href=REFERENCES.crown.url;a.textContent=REFERENCES.crown.title;a.target='_blank';a.rel='noopener noreferrer';get('inspect-sources').append(a);}
 highlight();updateLabels();
}
function evidenceSVG(){
 const unresolved=taxon!=='all'&&!observation(taxon,health).resolved;
 if(health==='healthy')return '<svg viewBox="0 0 600 245" role="img" aria-label="Healthy-site location unresolved; no spatial placement inferred"><path d="M52 45V199" stroke="#d6c49e" stroke-width="15"/><path d="M70 58H555V198H70Z" fill="#21373b"/><text x="300" y="115" text-anchor="middle" fill="#d8bc8b" font-size="22">Location unresolved</text><text x="300" y="145" text-anchor="middle" fill="#a9bcb4" font-size="13">No compatible healthy-site spatial evidence.</text><text x="300" y="180" text-anchor="middle" fill="#a9bcb4" font-size="12">Disease-specimen positions are not transferred.</text></svg>';
 let glyphs='';
 const draw=(id:string,x:number,y:number,k:number)=>{const t=TAXA.find(t=>t.id===id)!,opacity=taxon==='all'||taxon===id?1:.14;
  if(id==='ac')return `<path d="M${x-12} ${y+8}l17 -16m-10 6l-10 -8m14 6l10 8" stroke="${t.color}" stroke-width="4" stroke-linecap="round" opacity="${opacity}" fill="none"/>`;
  if(['pg','pi','pm'].includes(id))return `<g fill="${t.color}" opacity="${opacity}">${[[-8,-5],[1,-9],[8,-2],[-4,5],[6,7]].map(([a,b])=>`<ellipse cx="${x+a}" cy="${y+b}" rx="4" ry="${id==='pm'?4:5}"/>`).join('')}</g>`;
  return `<rect x="${x-16}" y="${y-3}" width="32" height="6" rx="3" fill="${t.color}" opacity="${opacity}" transform="rotate(${k%2?24:-22} ${x} ${y})"/>`;
 };
 for(let i=0;i<6;i++)glyphs+=draw('ac',100+(i%2)*35,85+Math.floor(i/2)*38,i);
 for(let i=0;i<6;i++)glyphs+=draw(i%2?'tf':'fn',230+(i%2)*65,86+Math.floor(i/2)*40,i);
 for(const [id,x,y,k]of [['pg',373,91,0],['pi',416,156,0],['pm',478,109,0]] as [string,number,number,number][])glyphs+=draw(id,x,y,k);
 if(unresolved)glyphs+='<rect x="168" y="96" width="309" height="58" rx="6" fill="#11262b"/><text x="322" y="120" text-anchor="middle" fill="#d8bc8b" font-size="17">Selected location unresolved</text><text x="322" y="141" text-anchor="middle" fill="#a9bcb4" font-size="11">Other source observations remain faintly visible.</text>';
 return `<svg viewBox="0 0 600 245" role="img" aria-label="Schematic across-biofilm observations from selected human periodontitis specimens; glyphs are not measured cells"><defs><linearGradient id="biofilm-matrix"><stop stop-color="#324740"/><stop offset="1" stop-color="#213b40"/></linearGradient></defs><path d="M52 50V198" stroke="#d6c49e" stroke-width="15"/><path d="M70 61Q184 42 300 65T532 62L535 191Q400 200 292 186T70 198Z" fill="url(#biofilm-matrix)"/>${glyphs}<text x="120" y="26" fill="#c8d8ce" text-anchor="middle" font-size="13">Tooth-facing</text><text x="269" y="26" fill="#c8d8ce" text-anchor="middle" font-size="13">Intermediate</text><text x="425" y="26" fill="#c8d8ce" text-anchor="middle" font-size="13">Microcolonies within biofilm</text><path d="M76 224H527m-8 -5l8 5l-8 5" stroke="#a2b9ab" fill="none"/><text x="300" y="218" text-anchor="middle" fill="#a2b9ab" font-size="11">Across attached plaque → fluid-side boundary</text></svg>`;
}
// Detail orientation: where the across-plaque view sits on the tooth. It marks the
// compartment only; species positions come from the evidence panel beside it.
function whereSVG(){
 const deep=health==='periodontitis',bottom=deep?96:70;
 return `<svg class="lens-where" viewBox="0 0 150 150" role="img" aria-label="Detail location: distal root surface of FDI 36 below the gum margin; tooth surface on the left of the detail, pocket fluid on the right"><path d="M38 40Q36 18 52 14Q66 22 74 14Q88 22 96 14Q114 18 110 40L108 56H40Z" fill="#efe7d6"/><path d="M40 56H108L106 70Q104 110 98 140Q92 140 90 112L84 78H64L58 112Q56 140 50 140Q44 110 42 70Z" fill="#dcc29a"/><path d="M106 ${deep?62:52}Q124 ${deep?60:50} 140 58V150H114Q106 120 106 ${deep?62:52}Z" fill="#d98278" opacity=".85"/><path d="M8 ${deep?60:50}Q24 ${deep?58:48} 40 52V150H8Z" fill="#d98278" opacity=".85"/><rect x="104" y="${deep?62:54}" width="7" height="${bottom-(deep?62:54)}" rx="2" fill="#d8bc8b" stroke="#fff4cc" stroke-width="1.2"/><text x="75" y="9" text-anchor="middle" fill="#a9bcb4" font-size="9">crown ↑</text><text x="128" y="${bottom+12}" text-anchor="middle" fill="#fff4cc" font-size="8">detail</text></svg>`;
}
function renderLens(){
 const visible=layer==='species'&&scale!=='mouth'&&scale!=='face'&&scale!=='crevice';lens.hidden=!visible;stage.classList.toggle('detail',scale==='biofilm');stage.classList.toggle('with-lens',visible&&scale!=='biofilm');
 get('detail').textContent=scale==='biofilm'?'← Return to pocket':'Explore biofilm detail →';
 const run=saved.cases.find(c=>c.id===health)!;lens.dataset.basis=basis;
 if(basis==='evidence'){
  get('lens-title').textContent='Across attached plaque';get('lens-note').textContent=health==='healthy'?'Healthy-site locations unresolved; disease-specimen observations are not transferred.':'Qualitative human observations · schematic glyphs, not counts or measured positions.';
  get('lens-content').innerHTML=whereSVG()+evidenceSVG();get('lens-axis').textContent='Across-biofilm organization ≠ bands down the pocket. No physical scale.';
 }else{
  get('lens-title').textContent='Assumed model coordinates';get('lens-note').textContent='Original 32 × 72 grid · fixed saved endpoint · no anatomical projection.';
  const c=document.createElement('canvas');c.width=run.nx;c.height=run.ny;c.setAttribute('aria-label','Saved '+health+' '+taxon+' density grid, original coordinates');c.setAttribute('role','img');c.dataset.case=run.id;c.dataset.taxon=taxon;c.dataset.sourceSha256=saved.source.sha256;c.dataset.frame=String(run.frame_index);
  const ctx=c.getContext('2d')!,pixels=ctx.createImageData(run.nx,run.ny);pixels.data.set(gridRGBA(run,taxon));ctx.putImageData(pixels,0,0);
  const image=document.createElement('img');image.id='grid-image';image.src=c.toDataURL('image/png');image.alt=c.getAttribute('aria-label')!;image.width=run.nx;image.height=run.ny;for(const [key,value]of Object.entries(c.dataset))image.dataset[key]=value;const badge=document.createElement('p');badge.className='lens-badge';badge.textContent='Saved teaching field · original model grid, not placed on the tooth';const layout=document.createElement('div');layout.className='model-detail';const chart=document.createElement('figure');chart.className='density-map';const entry=document.createElement('figcaption');entry.textContent='Model channel entrance · row 0';const end=document.createElement('figcaption');end.textContent='↓ Assumed axial depth · row 71';chart.append(entry,image,end);const forms=document.createElement('div');forms.innerHTML=morphologyPanel(taxon);layout.append(chart,forms);get('lens-content').replaceChildren(badge,layout);
  get('lens-axis').textContent='↓ Model axial depth · → Across model channel; not physical plaque thickness. Contrast normalized per selection.';
 }
}
function renderEnvironment(){
 const compartment=get<HTMLSelectElement>('compartment').value;
 const entries:Array<[string,string]>=[];
 if(compartment==='lumen'&&health==='periodontitis')entries.push(['Pocket-base oxygen · source cohort','Reported mean 13.3 mmHg, observed range 5–27 mmHg in untreated advanced human pockets. This is a source observation, not a value assigned to this illustrated site.']);
 else entries.push(['Oxygen at this location','Unknown. No compatible site-specific measurement is imported.']);
 entries.push(['pH at this location','Unknown. No field or fixed pH is assigned.']);
 entries.push(['Fluid flow / tissue perfusion',compartment==='lumen'?'Site-specific GCF flux unknown. Early plaque/gingivitis fluid observations do not calibrate this advanced pocket.':'Unknown. Pocket-fluid values are not transferred into this compartment.']);
 entries.push(['Cytokines and immune-cell contacts','Unknown. No response multiplier, contact inference or immune prediction is applied.']);
 const dl=get('environment-content');dl.replaceChildren();for(const [key,value]of entries){const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=value;dl.append(dt,dd);}
 inspect('Local environment','Measurements belong to a specific host, site, compartment and method.',entries[0][1],compartment==='lumen'&&health==='periodontitis'?SOURCES.oxygen.method+' '+SOURCES.oxygen.limit:'No compatible site-specific environmental values are available. Unknown is not zero.','Source observations · site values unresolved','oxygen');
}
function refresh(){
 if(fallback&&scale==='crevice')scale='pocket';
 if(scale==='crevice'&&stage.dataset.creviceSelection!==health+':'+taxon){geometry();if(creviceView==='wall')view(true);}
 document.body.classList.toggle('crevice-view',scale==='crevice');stage.classList.toggle('crevice-mode',scale==='crevice');get('crevice-controls').hidden=scale!=='crevice';get('crevice-key').hidden=scale!=='crevice';
 for(const b of document.querySelectorAll<HTMLButtonElement>('[data-crevice-view]'))b.setAttribute('aria-pressed',String(b.dataset.creviceView===creviceView));
 stage.dataset.health=health;stage.dataset.scale=scale;stage.dataset.layer=layer;stage.dataset.basis=basis;stage.dataset.taxon=taxon;stage.dataset.cutaway=String(cutaway);stage.dataset.host='human';stage.dataset.anatomy=scale==='face'||scale==='mouth'?'BodyParts3D reference atlas; authored jaw opening; independent pocket detail':'Reference-sized crown; surrounding geometry authored; expert review pending';stage.dataset.referenceScale=JSON.stringify(CROWN_REFERENCE);stage.dataset.neuro=String(neuro);
 get('neuro').setAttribute('aria-pressed',String(neuro));get<HTMLButtonElement>('neuro').disabled=scale==='crevice';
 if(!(neuro&&['nerve','vessels'].includes(selected)))get('nerve-inspector').hidden=true;
 for(const [attr,value]of [['health',health],['scale',scale==='biofilm'?'pocket':scale],['layer',layer],['basis',basis]])for(const b of document.querySelectorAll<HTMLButtonElement>(`[data-${attr}]`))b.setAttribute('aria-pressed',String(b.dataset[attr]===value));
 get('labels').setAttribute('aria-pressed',String(showLabels));get('cutaway').setAttribute('aria-pressed',String(cutaway));get<HTMLButtonElement>('cutaway').disabled=scale==='mouth'||scale==='face'||scale==='biofilm'||scale==='crevice'||fallback;
 for(const b of document.querySelectorAll<HTMLButtonElement>('[data-scale="mouth"],[data-scale="face"],[data-scale="crevice"]'))b.disabled=fallback;
 const copy=sceneCopy(scale,health,cutaway);get('scene-title').textContent=copy.title;get('scene-subtitle').textContent=copy.subtitle;
 if(scale==='pocket'||scale==='tooth')get('scene-subtitle').textContent='Reference-sized crown · pocket and tissue dimensions authored · '+(neuro?'neurovascular diameters enlarged':'no calibrated pocket ruler');
 get('orientation').hidden=scale==='biofilm';get('structure-controls').hidden=layer!=='structure';get('species-controls').hidden=layer!=='species';document.querySelector<HTMLElement>('.basis-controls')!.hidden=layer!=='species';get('environment-controls').hidden=layer!=='environment';get('environment-content').hidden=layer!=='environment';
 get<HTMLSelectElement>('taxon').value=taxon;
 for(const b of get('legend').querySelectorAll<HTMLButtonElement>('button'))b.setAttribute('aria-pressed',String(b.dataset.taxon===taxon));
 const note=get('scene-note');
 if(layer==='species'){
  note.textContent=basis==='model'?'Assumed model coordinates · Dominant color identifies the largest modeled density in a bin; other taxa may overlap. Saved fields have unresolved numerical and biological limitations.':health==='healthy'?'Location unresolved for healthy anatomy. The selected spatial evidence concerns advanced periodontitis; no disease-specimen positions are copied.':'Source-scoped observations across attached biofilm. These schematic arrangements do not establish universal pocket-depth positions, cell contacts or abundance.';
  inspectSpecies();
 }else if(layer==='environment'){note.textContent='Source observations remain separate from site-specific values. No continuous oxygen/pH field is invented.';renderEnvironment();}
 else {note.textContent='Margin, epithelial attachment and bone crest are separate authored presets. The crown is reference-sized; pocket and tissue geometry remain authored. Expert anatomical review is pending.';inspect(names[selected],...descriptions[selected],SOURCES.anatomy.method+' '+SOURCES.anatomy.limit,'Authored anatomy · expert review pending','anatomy');}
 if(scale==='mouth'||scale==='face'){get('scene-title').textContent=scale==='face'?'Craniofacial atlas':'Oral atlas · opened display pose';get('scene-subtitle').textContent='BodyParts3D adult male reference · pocket detail is an independent authored specimen';get('orientation').textContent='Subject left = +x · crown ↑';note.textContent=ATLAS_CREDIT+'. Jaw opening is authored; no disease geometry or species map is inferred for this atlas.';inspectAtlas(display.find(m=>m.userData.atlasIdentity?.fdi===36)!.userData.atlasIdentity);}else get('orientation').textContent='Crown ↑ / apex ↓';
 if(scale==='crevice'){
  get('scene-title').textContent='Inside the tooth–gum crevice';get('scene-subtitle').textContent='Magnified schematic · '+(health==='healthy'?'Healthy sulcus':'Periodontitis pocket')+' · width and cell size exaggerated';get('orientation').textContent=creviceView==='wall'?'Tooth-facing surface':'Entrance ↑ · attachment ↓';get('crevice-key').querySelector('strong')!.textContent=creviceView==='wall'?'Plaque surface · margin and support retained':'Tooth surface · fluid space · tissue lining';stage.dataset.creviceView=creviceView;note.textContent=CREVICE_CLAIM;
  if(layer==='species'){
   if(taxon==='all'){inspect('Nine microbial forms','An illustrative 3D form library beside the tooth-facing plaque.','All nine exemplar identities remain in the scene. Positions, counts and relative sizes are display choices, not colonization or measured cell contacts. Click a cell or select a taxon to inspect its shape.','Source-linked qualitative microscopy. Healthy/disease presets do not infer which organisms exist at this site.','3D forms · placement unresolved','plaque');get('inspect-sources').replaceChildren();const a=document.createElement('a');a.href=SOURCES.plaque.url;a.textContent='Attached-plaque compartment reference ↗';a.target='_blank';a.rel='noopener noreferrer';get('inspect-sources').append(a);}
   else {const t=TAXA.find(t=>t.id===taxon)!,m=MORPHOLOGY[taxon];inspect(t.name,m.observation,'A selectable 3D '+m.label.toLowerCase()+'. Other exemplar taxa remain in the scene in muted colors. '+CREVICE_CLAIM,m.context,'Qualitative morphology · location unresolved','plaque');get('inspect-sources').replaceChildren();const a=document.createElement('a');a.href=m.source;a.textContent='Morphology microscopy reference ↗';a.target='_blank';a.rel='noopener noreferrer';get('inspect-sources').append(a);}
  }else if(layer==='structure')get('inspect-shown').textContent=descriptions[selected][1]+' This magnified crevice is a separately authored explanatory surface; no measured registration or millimetre scale is implied.';
 }
 renderLens();if(!get('storyboard').hidden)renderStory();
 if(neuro&&['nerve','vessels'].includes(selected))renderNeuralInspector();resize();highlight();
}
function setScale(next:Scale){
 if(next!=='biofilm'&&!get('storyboard').hidden)closeStory();
 stopCameraTravel();remember();const previous=scale;scale=next;
 if(next!=='biofilm'){geometry();view();}else if(previous==='mouth'||previous==='face'){geometry();view();}
 refresh();resize();
}
function resize(){if(renderer){const h=scale==='crevice'?stage.clientHeight-175:stage.classList.contains('with-lens')?stage.clientHeight-247:stage.clientHeight-62;renderer.setPixelRatio(Math.min(devicePixelRatio,1.5));stage.dataset.pixelRatio=String(renderer.getPixelRatio());renderer.setSize(stage.clientWidth,h,false);camera.aspect=stage.clientWidth/h;camera.updateProjectionMatrix();composer?.setPixelRatio(renderer.getPixelRatio());composer?.setSize(stage.clientWidth,h);}updateLabels();dirty=true;}
for(const b of document.querySelectorAll<HTMLButtonElement>('[data-health]'))b.onclick=()=>{health=b.dataset.health as Health;geometry();if(scale==='crevice')view(true);refresh();};
for(const b of document.querySelectorAll<HTMLButtonElement>('[data-scale]'))b.onclick=()=>setScale(b.dataset.scale as Scale);
for(const b of document.querySelectorAll<HTMLButtonElement>('[data-layer]'))b.onclick=()=>{layer=b.dataset.layer as Layer;if(layer!=='species'&&scale==='biofilm')setScale('pocket');refresh();};
for(const b of document.querySelectorAll<HTMLButtonElement>('[data-basis]'))b.onclick=()=>{closeStory();basis=b.dataset.basis as Basis;if(basis==='model'&&scale!=='biofilm')setScale('biofilm');else refresh();};
for(const b of document.querySelectorAll<HTMLButtonElement>('[data-crevice-view]'))b.onclick=()=>{creviceView=b.dataset.creviceView as CreviceView;remembered.delete('crevice');refresh();travelCreviceCamera();};
get('lens-content').addEventListener('click',e=>{const b=(e.target as HTMLElement).closest<HTMLButtonElement>('[data-morphology]');if(b){taxon=b.dataset.morphology!;refresh();}});
const select=get<HTMLSelectElement>('taxon');
const tissueSelect=get<HTMLSelectElement>('tissue');for(const id of ['enamel','dentin','pulp','cementum','pdl','gingiva','bone','supragingival','plaque','lumen','epithelium','connective','vessels','nerve','cheek','tongue'] as PocketId[]){const option=document.createElement('option');option.value=id;option.textContent=names[id];tissueSelect.append(option);}tissueSelect.value='plaque';tissueSelect.onchange=()=>selectTissue(tissueSelect.value as PocketId);
function renderNeuralInspector(){
 const index=Number(get<HTMLSelectElement>('nerve-step').value)||0;
 get('nerve-path').replaceChildren();
 NERVE_STEPS.forEach((step,i)=>{const li=document.createElement('li');li.textContent=step.title;li.dataset.active=String(i===index);get('nerve-path').append(li);});
 get('nerve-description').textContent=NERVE_STEPS[index].text;
 get('inspect-title').textContent=selected==='vessels'?'Pulpal neurovascular supply':'Tooth → trigeminal sensory pathway';
 get('inspect-is').textContent='Named anatomical connections, with pulpal and gingival distributions kept separate.';
 get('inspect-shown').textContent=scale==='mouth'||scale==='face'?'This atlas subset has no identified inferior alveolar nerve or pulpal vessels. The named pathway is explained here; open Tooth to see a separately authored schematic route.':'Gold nerve, red artery and blue vein. Coordinates and diameters are schematic; central connections stop at brainstem entry.';
 get('inspect-support').textContent=REFERENCES.nerves.population+'. '+REFERENCES.nerves.limit;
 get('inspect-tag').textContent='Anatomical pathway · rendered course illustrative';
 get('inspect-sources').replaceChildren();
 for(const ref of [REFERENCES.nerves,REFERENCES.pathway]){const a=document.createElement('a');a.href=ref.url;a.textContent=ref.title+' ↗';a.target='_blank';a.rel='noopener noreferrer';get('inspect-sources').append(a);}
}
for(const [i,step]of NERVE_STEPS.entries()){const o=document.createElement('option');o.value=String(i);o.textContent=`${i+1}. ${step.title}`;get('nerve-step').append(o);}
get<HTMLSelectElement>('nerve-step').onchange=renderNeuralInspector;
get('neuro').onclick=()=>{neuro=!neuro;get('neuro').setAttribute('aria-pressed',String(neuro));geometry();if(neuro){selected='nerve';showInspector();get('nerve-inspector').hidden=false;renderNeuralInspector();}else get('nerve-inspector').hidden=true;refresh();dirty=true;};
for(const p of PARAMETERS){
 const tr=document.createElement('tr'),ref=REFERENCES[p.source];
 for(const text of [p.name+' — '+p.definition,p.reference,p.chosen+' · '+p.units]){const td=document.createElement('td');td.textContent=text;tr.append(td);}
 const td=document.createElement('td'),tag=document.createElement('strong'),a=document.createElement('a'),detail=document.createElement('p');tag.textContent=p.classification;a.href=ref.url;a.textContent=ref.title;a.target='_blank';a.rel='noopener noreferrer';detail.textContent=ref.population+'. '+ref.method+'. '+ref.limit;td.append(tag,a,detail);tr.append(td);get('reference-table').querySelector('tbody')!.append(tr);
}
for(const [name,description,key]of [
 ['Supragingival plaque','Tooth surface above the gingival margin. Organized consortia have been imaged here; those arrangements are not transferred into a deep pocket.','supra'],
 ['Subgingival plaque','Tooth-associated plaque below the margin. Selected disease-specimen observations are shown separately from healthy sites and saved model grids.','plaque'],
 ['Tongue dorsum','A shedding epithelial habitat with structured consortia. Static images support inferred dynamics, not a measured growth film.','tongue'],
 ['Mucosal habitats and saliva','Oral sites have distinct community compositions. Saliva connects habitats but its composition is not a resident-site position map.','habitats']
 ] as const){const article=document.createElement('article'),h=document.createElement('h3'),p=document.createElement('p'),a=document.createElement('a');h.textContent=name;p.textContent=description;a.href=REFERENCES[key].url;a.textContent=REFERENCES[key].title;a.target='_blank';a.rel='noopener noreferrer';article.append(h,p,a);get('habitat-records').append(article);}
// The storyboard lives in the same enlarged detail scene. It never writes to
// anatomy presets, saved grids or scientific records.
lens.insertBefore(get('storyboard'),get('detail'));
function stopStory(){playing=false;storyLast=0;get('story-play').setAttribute('aria-pressed','false');get('story-play').textContent='Play stages';}
function renderStory(){
 const position=reducedMotion?Math.floor(storyProgress):storyProgress;
 const step=stageAt(sequence,position),ref=sourceForStage(sequence,position);
 get('lens-title').textContent='Biofilm story · illustrative stages';get('lens-note').textContent='Stage sequence, not a calibrated simulation. Species observations and saved fields remain separate.';
 get('story-art').innerHTML=sequenceSVG(sequence,position);get('story-title').textContent=step.title;get('story-status').textContent=step.status+' · animation illustrative';get('story-text').textContent=step.text;
 const link=get<HTMLAnchorElement>('story-source');link.href=ref.url;link.textContent=ref.title+' ↗';
 get<HTMLInputElement>('story-progress').value=String(storyProgress);
 stage.dataset.sequence=sequence;stage.dataset.storyStage=String(step.index);stage.dataset.storyProgress=String(storyProgress);
}
function closeStory(){stopStory();get('storyboard').hidden=true;get('lens-content').hidden=false;get('lens-axis').hidden=false;renderLens();}
get('open-story').onclick=()=>{layer='species';setScale('biofilm');get('storyboard').hidden=false;get('lens-content').hidden=true;get('lens-axis').hidden=true;get('lens-title').textContent='Biofilm story · illustrative stages';get('lens-note').textContent='Stage sequence, not a calibrated simulation. Species observations and saved fields remain separate.';renderStory();};
get('story-close').onclick=closeStory;
get('story-play').onclick=()=>{if(playing){stopStory();return;}if(storyProgress>=5)storyProgress=0;playing=true;storyLast=0;get('story-play').setAttribute('aria-pressed','true');get('story-play').textContent='Pause';};
get('story-reset').onclick=()=>{stopStory();storyProgress=0;renderStory();};
get<HTMLInputElement>('story-progress').oninput=()=>{stopStory();storyProgress=Number(get<HTMLInputElement>('story-progress').value);renderStory();};
get<HTMLSelectElement>('sequence').onchange=()=>{stopStory();sequence=get<HTMLSelectElement>('sequence').value as Sequence;storyProgress=0;renderStory();};
document.addEventListener('visibilitychange',()=>{if(document.hidden)stopStory();});
for(const t of TAXA){const option=document.createElement('option');option.value=t.id;option.textContent=t.name;select.append(option);const b=document.createElement('button'),swatch=document.createElement('i');swatch.style.background=t.color;b.dataset.taxon=t.id;b.setAttribute('aria-pressed','false');b.append(swatch,document.createTextNode(t.short));b.onclick=()=>{taxon=t.id;showInspector();refresh();};get('legend').append(b);}
select.onchange=()=>{taxon=select.value;showInspector();refresh();};
get('detail').onclick=()=>{closeStory();setScale(scale==='biofilm'?'pocket':'biofilm');};
get('labels').onclick=()=>{showLabels=!showLabels;refresh();};get('cutaway').onclick=()=>{cutaway=!cutaway;geometry();refresh();};get('reset').onclick=()=>{if(scale==='biofilm')setScale('pocket');remembered.delete(scale);view(true);};
get('close-inspector').onclick=()=>{inspector=false;get('inspector').hidden=true;get('workspace').classList.add('inspector-closed');get('show-inspector').hidden=false;resize();};get('show-inspector').onclick=showInspector;
get<HTMLSelectElement>('compartment').onchange=renderEnvironment;
get('provenance').textContent='Source: '+saved.source.file+' · SHA-256 '+saved.source.sha256+'. '+saved.source.precision+' Both cases show their fixed endpoint; model time is not physical elapsed time.';
get('render-check').onclick=()=>{if(!renderer){get('performance').textContent='2D fallback active; WebGL performance not measured.';return;}if(scale==='biofilm')setScale('pocket');measuring=true;measureStart=performance.now();frameCount=0;renderMS=[];get('performance').textContent='Measuring 60 rendered frames…';dirty=true;};
canvas.addEventListener('keydown',e=>{if(!controls)return;if(e.key==='Home'){e.preventDefault();view(true);return;}if(!e.key.startsWith('Arrow'))return;e.preventDefault();stopCameraTravel();const offset=camera.position.clone().sub(controls.target),s=new T.Spherical().setFromVector3(offset);s.theta+=e.key==='ArrowLeft'?-.10:e.key==='ArrowRight'?.10:0;s.phi=T.MathUtils.clamp(s.phi+(e.key==='ArrowUp'?-.08:e.key==='ArrowDown'?.08:0),.1,Math.PI-.1);camera.position.copy(controls.target).add(new T.Vector3().setFromSpherical(s));controls.update();dirty=true;updateLabels();});
let down=[0,0];canvas.addEventListener('pointerdown',e=>{stopCameraTravel();down=[e.clientX,e.clientY];});canvas.addEventListener('wheel',stopCameraTravel,{passive:true});canvas.addEventListener('pointerup',e=>{if(!renderer||Math.hypot(e.clientX-down[0],e.clientY-down[1])>5)return;const r=canvas.getBoundingClientRect(),pointer=new T.Vector2((e.clientX-r.left)/r.width*2-1,-(e.clientY-r.top)/r.height*2+1),ray=new T.Raycaster();ray.setFromCamera(pointer,camera);const hits=ray.intersectObjects(display).filter(h=>h.object.visible),hit=hits.find(h=>!h.object.userData.context)||hits[0];if(hit){if(hit.object.userData.taxon){taxon=hit.object.userData.taxon;layer='species';showInspector();refresh();}else if(hit.object.userData.atlasIdentity)inspectAtlas(hit.object.userData.atlasIdentity,true);else selectTissue(hit.object.userData.tissue)};});
function fallbackScene(reason:string){fallback=true;stage.dataset.renderer='2d-fallback';canvas.hidden=true;get('fallback').hidden=false;status(reason+' The 2D anatomy and all evidence/model controls remain available.');updateFallback();}
try{
 if(new URLSearchParams(location.search).get('render')==='2d')throw Error('2D fallback requested');
 renderer=new T.WebGLRenderer({canvas,antialias:true,alpha:true,powerPreference:'high-performance'});renderer.localClippingEnabled=true;renderer.setClearColor('#0c171a',0);renderer.toneMapping=T.NeutralToneMapping;renderer.toneMappingExposure=.74;renderer.shadowMap.enabled=true;renderer.shadowMap.type=T.PCFShadowMap;
 mats=materials();boneCapMap=boneSectionTexture();fiberCapMap=connectiveSectionTexture();for(const m of Object.values(mats))m.envMapIntensity=.28;mats.enamel.envMapIntensity=.42;scene.environment=environment(renderer);
 // Restrained three-point lighting: warm soft key, cool low fill, neutral rim for silhouette separation.
 const key=new T.DirectionalLight('#ffe9d6',2.1);key.position.set(-15,24,18);key.castShadow=true;key.shadow.mapSize.set(2048,2048);key.shadow.camera.left=-25;key.shadow.camera.right=25;key.shadow.camera.top=25;key.shadow.camera.bottom=-25;key.shadow.bias=-.0002;key.shadow.normalBias=.02;key.shadow.radius=3;
 const fill=new T.DirectionalLight('#b4cfe0',.55);fill.position.set(18,2,16);const rim=new T.DirectionalLight('#e8eef2',1.25);rim.position.set(-8,12,-18);scene.add(key,fill,rim,new T.HemisphereLight('#e4e2dc','#1a2124',.22));
 controls=new OrbitControls(camera,canvas);controls.enableDamping=false;controls.enablePan=false;controls.minDistance=20;controls.maxDistance=600;controls.addEventListener('change',()=>{dirty=true;updateLabels();});
 composer=new EffectComposer(renderer);composer.addPass(new RenderPass(scene,camera));
 ao=new SSAOPass(scene,camera,512,512,16);ao.kernelRadius=.8;ao.minDistance=.0002;ao.maxDistance=.025;composer.addPass(ao);composer.addPass(new OutputPass());
 stage.dataset.renderEffects='Contact shadows · physical materials · filmic tone mapping';
 const gl=renderer.getContext(),gpu=gl.getExtension('WEBGL_debug_renderer_info');stage.dataset.gpuRenderer=String(gl.getParameter(gpu?gpu.UNMASKED_RENDERER_WEBGL:gl.RENDERER));
 stage.dataset.renderer='webgl2';geometry();resize();view(true);status('Local 3D scene ready · source-scoped evidence · expert anatomical review pending.');
 canvas.addEventListener('webglcontextlost',e=>{e.preventDefault();fallbackScene('WebGL context was lost.');refresh();});
}catch(error){fallbackScene((error as Error).message+'.');}
refresh();resize();if(controls)view(true);document.body.dataset.ready='true';
new ResizeObserver(resize).observe(stage);
function tick(){requestAnimationFrame(tick);if(cameraTravel&&controls){const t=T.MathUtils.clamp((performance.now()-cameraTravel.start)/460,0,1),ease=t*t*(3-2*t);camera.position.lerpVectors(cameraTravel.from,cameraTravel.to,ease);controls.target.lerpVectors(cameraTravel.fromTarget,cameraTravel.toTarget,ease);controls.update();dirty=true;if(t===1)stopCameraTravel();}if(playing){const now=performance.now();if(!storyLast)storyLast=now;storyProgress=Math.min(5,storyProgress+Math.min(100,now-storyLast)/2600);storyLast=now;renderStory();if(storyProgress>=5)stopStory();}
 if(renderer&&!fallback&&(dirty||measuring)&&scale!=='biofilm'){
 const start=performance.now();composer?composer.render():renderer.render(scene,camera);dirty=false;updateLabels();
 if(measuring){renderMS.push(performance.now()-start);frameCount++;if(frameCount>=60){measuring=false;const wall=performance.now()-measureStart,sorted=renderMS.slice().sort((a,b)=>a-b),fps=frameCount*1000/wall;stage.dataset.measuredFps=fps.toFixed(1);stage.dataset.renderMedianMs=sorted[Math.floor(sorted.length/2)].toFixed(2);get('performance').textContent=` ${fps.toFixed(1)} frames/s over ${frameCount} frames · median CPU submission ${sorted[Math.floor(sorted.length/2)].toFixed(2)} ms · ${Math.round(stage.clientWidth)} × ${Math.round(stage.clientHeight)} CSS pixels. This measures display execution, not biological accuracy.`;}}
 }}tick();
