import * as T from 'three';
import {EffectComposer} from 'three/addons/postprocessing/EffectComposer.js';
import {RenderPass} from 'three/addons/postprocessing/RenderPass.js';
import {SSAOPass} from 'three/addons/postprocessing/SSAOPass.js';
import {OutputPass} from 'three/addons/postprocessing/OutputPass.js';
import {OrbitControls} from 'three/addons/controls/OrbitControls.js';
import {toCreasedNormals} from 'three/addons/utils/BufferGeometryUtils.js';
import {environment} from './materials';
import {pocketMaterials as materials,boneSectionTexture,connectiveSectionTexture} from './pocket-materials';
import {quarterSection} from './pocket-section';
import {buildPocket,buildOralContext,buildOralSides,ANCHORS,PROFILES,type PocketId,type PocketPart} from './pocket-geometry';
import {TAXA,SOURCES,PRESETS,observation,validateSaved,gridRGBA,type Health,type Scale,type Layer,type Basis,type SavedInput} from './pocket-science';
import rawInput from './pocket-input.json';

const get=<E extends HTMLElement=HTMLElement>(id:string)=>document.getElementById(id) as E;
const stage=get('stage'),canvas=get<HTMLCanvasElement>('scene'),labels=get('label-layer'),lens=get('lens');
const saved=validateSaved(rawInput as SavedInput);
const entryView=new URLSearchParams(location.search);
let health:Health=entryView.get('health')==='healthy'?'healthy':'periodontitis',scale:Scale='pocket',layer:Layer='species',basis:Basis='evidence',taxon='all',selected:PocketId='plaque',cutaway=entryView.get('cutaway')!=='0',showLabels=true,inspector=true;
let renderer:T.WebGLRenderer|undefined,controls:OrbitControls|undefined,mats:ReturnType<typeof materials>|undefined,dirty=true,fallback=false;
const scene=new T.Scene(),camera=new T.PerspectiveCamera(34,1,.1,350),root=new T.Group();scene.add(root);
const connector=document.createElementNS('http://www.w3.org/2000/svg','svg');connector.classList.add('model-connector');connector.setAttribute('aria-hidden','true');stage.append(connector);
let boneCapMap:T.Texture|undefined,fiberCapMap:T.Texture|undefined;
let composer:EffectComposer|undefined,ao:SSAOPass|undefined;
const cache=new Map<string,PocketPart[]>(),display:T.Mesh[]=[],capMats=new Map<string,T.MeshPhysicalMaterial>();
const tissueAnchors=new Map<PocketId,T.Vector3>();
const remembered=new Map<Scale,{position:number[];target:number[]}>();
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
 lumen:['The fluid space between the tooth-facing surface and the epithelial wall.','A slit-like opening and an enlarged illustrative interior. Its width is unknown; it is not filled with stretched microbial layers.'],
 epithelium:['The tooth-associated epithelial boundary; the junctional attachment is distinct from a clinical probing endpoint.','Sulcular lining and an apical attachment region are shown as a separate tissue. Boundaries are not measured histology.'],
 connective:['Vascular connective tissue beneath the epithelial boundary.','A distinct cut-face region behind the lining. No immune migration, inflammation dynamics or tissue oxygen field is computed.']
};
function status(text:string){get('status').textContent=text;}
function remember(){if(controls&&scale!=='biofilm')remembered.set(scale,{position:camera.position.toArray(),target:controls.target.toArray()});}
function view(reset=false){
 if(!controls)return;
 const stored=!reset&&remembered.get(scale);if(stored){camera.position.fromArray(stored.position);controls.target.fromArray(stored.target);}
 else if(scale==='mouth'){camera.position.set(38,34,142);controls.target.set(0,12,10);}
 else{
  const box=modelBox();if(scale==='pocket')box.min.y=-10;const center=box.getCenter(new T.Vector3());
  const direction=new T.Vector3(.48,.25,1).normalize(),right=new T.Vector3().crossVectors(camera.up,direction).normalize(),up=new T.Vector3().crossVectors(direction,right).normalize();
  const tangent=Math.tan(T.MathUtils.degToRad(camera.fov/2)),fraction=stage.clientWidth<450?.43:.66;let distance=20;
  const point=new T.Vector3();
  for(const mesh of display){if(mesh.userData.context)continue;const p=mesh.geometry.getAttribute('position');for(let i=0;i<p.count;i++){
   point.fromBufferAttribute(p,i);if(scale==='pocket'&&point.y<-10)continue;point.sub(center);distance=Math.max(distance,Math.abs(point.dot(up))/(tangent*.92)+point.dot(direction),Math.abs(point.dot(right))/(tangent*camera.aspect*fraction)+point.dot(direction));
  }}
  controls.target.copy(center);camera.position.copy(center).addScaledVector(direction,distance*(scale==='tooth'?1.10:1));
 }
 controls.update();dirty=true;updateLabels();
}
// Translucent oral-side context never drives framing or label placement.
function modelBox(){const box=new T.Box3();for(const mesh of display)if(!mesh.userData.context)box.expandByObject(mesh);return box;}
function material(id:PocketId){
 if(id in mats!)return mats![id as keyof typeof mats];
 const colors={plaque:'#b9a46a',supragingival:'#d4c79a',lumen:'#2c5862',epithelium:'#e9a397',connective:'#dba196',cheek:'#c98a8a'};
 if(id==='cheek')return new T.MeshPhysicalMaterial({color:colors.cheek,roughness:.5,clearcoat:.2,transparent:true,opacity:.22,depthWrite:false,side:T.DoubleSide});
 return new T.MeshPhysicalMaterial({color:colors[id as keyof typeof colors],roughness:id==='lumen'?.3:.6,clearcoat:id==='lumen'?.3:.1,side:T.DoubleSide});
}
const extraMats=new Map<PocketId,T.MeshPhysicalMaterial>();
let tongueContext:T.MeshPhysicalMaterial|undefined;
function tongueMat(){return tongueContext??=new T.MeshPhysicalMaterial({color:'#c4777b',roughness:.45,clearcoat:.3,transparent:true,opacity:.5,depthWrite:false,side:T.DoubleSide});}
function mat(id:PocketId){if(id in mats!)return material(id);if(!extraMats.has(id))extraMats.set(id,material(id));return extraMats.get(id)!;}
function clear(){for(const mesh of display){root.remove(mesh);mesh.geometry.dispose();}display.length=0;}
function anchorTarget(id:PocketId):[number,number,number]{
 if(id==='bone'&&!cutaway)return [-6.4,-8.0,4.6];
 const anchor=[...ANCHORS[id]!] as [number,number,number];if(['plaque','lumen','epithelium'].includes(id)){
   anchor[1]=id==='epithelium'?PRESETS[health].attachment+.10:(PRESETS[health].margin+PRESETS[health].attachment)/2;
   const y=anchor[1],t=Math.max(0,Math.min(1,(-y-2)/11.60));
   const neck=y>=-4.4?4.65*Math.sqrt(Math.max(0,1-((y+1.2)/3.2)**2)):0;
   const rootX=y<=-2?2+1.05*Math.sin(t*Math.PI/2)+.85*t*t+1.74*(1-.5*t)*(1-t**3.2)**.42+.03:0;
   anchor[0]=Math.max(neck,rootX)+.13+(id==='plaque'?.035:id==='lumen'?.24:.44);anchor[2]=0;
 }
 return anchor;
}
function geometry(){
 updateFallback();if(!renderer||!mats)return;
 clear();const key=scale==='mouth'?'mouth':health;
 // Split normals at Boolean contact creases so shading does not smear across tissue edges.
 if(!cache.has(key))cache.set(key,scale==='mouth'?buildOralContext():buildPocket(health).map(p=>{const g=toCreasedNormals(p.geometry,T.MathUtils.degToRad(48));p.geometry.dispose();return {...p,geometry:g};}));
 let open=0,loops=0,triangles=0;
 if(scale==='tooth'&&!cache.has('sides'))cache.set('sides',buildOralSides());
 const context=scale==='tooth'?cache.get('sides')!:[];
 for(const part of [...cache.get(key)!,...context]){const isContext=context.includes(part);
  let g:T.BufferGeometry,cap:T.BufferGeometry|undefined;
  if(cutaway&&scale!=='mouth'){const clipped=quarterSection(part.geometry);g=clipped.surface;cap=clipped.cap;if(!isContext){open+=clipped.openChains;loops+=clipped.loops;}}
  else g=part.geometry.clone();
  if(part.id==='enamel'){
   // Optical cervical warmth, authored art rather than a measured mineral map.
   for(const mesh of [g,cap].filter(Boolean) as T.BufferGeometry[]){const p=mesh.getAttribute('position'),colors=[];for(let i=0;i<p.count;i++){const tint=new T.Color('#ddd0b4').lerp(new T.Color('#f2eee6'),T.MathUtils.clamp((p.getY(i)+.3)/5.2,0,1)**.8);colors.push(tint.r,tint.g,tint.b);}mesh.setAttribute('color',new T.Float32BufferAttribute(colors,3));}
   mat('enamel').vertexColors=true;
  }
  const mesh=new T.Mesh(g,isContext&&part.id==='tongue'?tongueMat():mat(part.id));mesh.name=part.name;mesh.userData={tissue:part.id,context:isContext};mesh.castShadow=!isContext;mesh.receiveShadow=true;if(isContext)mesh.renderOrder=30;root.add(mesh);display.push(mesh);triangles+=(g.index?.count||g.getAttribute('position').count)/3;
  if(cap){if(cap.getAttribute('position').count){
   const capKey=part.id==='bone'?part.name:part.id;
   if(!capMats.has(capKey)){
    const m=mat(part.id).clone();
    m.map=part.id==='bone'?(part.name.includes('Trabecular')?boneCapMap!:null):part.id==='connective'?fiberCapMap!:mat(part.id).map;
    if(['gingiva','connective'].includes(part.id))m.color.set('#ffffff');
    m.bumpMap=part.id==='bone'?m.map:null;m.bumpScale=.04;m.clearcoat=.04;m.roughness=.76;
    if(part.id==='dentin')m.color.set('#f3e2c6');if(part.id==='bone')m.color.set(part.name.includes('Trabecular')?'#ffffff':'#d8cfbd');if(part.id==='pulp')m.color.set('#f0d4cc');if(part.id==='enamel')m.color.set('#f4eee2');m.depthWrite=true;m.envMapIntensity=.14;m.polygonOffset=true;m.polygonOffsetFactor=-rank[part.id];m.polygonOffsetUnits=-rank[part.id];m.userData.tissue=part.id;capMats.set(capKey,m);
   }
   const surface=new T.Mesh(cap,capMats.get(capKey));surface.name=part.name+' section';surface.userData={tissue:part.id};surface.renderOrder=10+rank[part.id];root.add(surface);display.push(surface);triangles+=cap.getAttribute('position').count/3;
  }else cap.dispose();}
 }
 if(scale==='mouth'){
  const ring=new T.Mesh(new T.TorusGeometry(3.8,.18,12,64),new T.MeshBasicMaterial({color:'#c8dcba'}));ring.rotation.x=-Math.PI/2;ring.position.set(24.4,7,9.2);ring.userData={tissue:'plaque',focus:true};root.add(ring);display.push(ring);
 }
 tissueAnchors.clear();const anchorCosts=new Map<PocketId,number>();
 for(const mesh of display){const id=mesh.userData.tissue as PocketId;if(!ANCHORS[id])continue;const target=new T.Vector3(...anchorTarget(id)),p=mesh.geometry.getAttribute('position'),point=new T.Vector3();
  for(let i=0;i<p.count;i++){point.fromBufferAttribute(p,i);const d=point.distanceToSquared(target);if(d<(anchorCosts.get(id)??Infinity)){anchorCosts.set(id,d);tissueAnchors.set(id,point.clone());}}
 }
 stage.dataset.tissueAnchors=JSON.stringify(Object.fromEntries([...tissueAnchors].map(([id,p])=>[id,p.toArray()])));
 stage.dataset.openContours=String(open);stage.dataset.closedContours=String(loops);stage.dataset.triangles=String(Math.round(triangles));
 highlight();dirty=true;updateLabels();
 if(open)status('Some section contours are open; geometry requires correction.');
}
function highlight(){
 if(!mats)return;
 for(const id of Object.keys(rank) as PocketId[]){const m=mat(id);m.emissive.set(id===selected?'#231407':'#000000');m.emissiveIntensity=.08;for(const cap of capMats.values())if(cap.userData.tissue===id){cap.emissive.copy(m.emissive);cap.emissiveIntensity=.06;}}
 if(layer==='species'){mat('plaque').emissive.set('#a17c36');mat('plaque').emissiveIntensity=.28;}
 dirty=true;
}
function updateFallback(){
 const v=PRESETS[health];get('fallback-lumen').setAttribute('d',PROFILES[health].lumen);get('fallback-gum').setAttribute('d',PROFILES[health].gingiva);
 get('fallback-bone').setAttribute('d',`M10.2,${v.crest} L5.05,${v.crest} Q5.6,-11 5,-14.5 L2,-16.2 L10.3,-16.2 Z`);
 get('fallback-plaque').setAttribute('d',`M4.45,${v.margin-.2} L4.55,${v.margin-.2} L4.6,${v.attachment} L4.48,${v.attachment} Z`);
}
const labelIds:PocketId[]=['enamel','dentin','pulp','gingiva','supragingival','cementum','plaque','lumen','epithelium','pdl','bone','cheek','tongue'];
const labelElements=new Map<PocketId,HTMLButtonElement>();
const leaders=document.createElementNS('http://www.w3.org/2000/svg','svg');leaders.setAttribute('class','leaders');leaders.setAttribute('aria-hidden','true');labels.append(leaders);
for(const id of labelIds){const b=document.createElement('button');b.textContent=names[id];b.setAttribute('aria-label','Inspect '+names[id]);b.addEventListener('click',()=>selectTissue(id));labels.append(b);labelElements.set(id,b);}
function updateLabels(){
 camera.updateMatrixWorld();
 if(controls){stage.dataset.camera=JSON.stringify(camera.position.toArray());stage.dataset.cameraTarget=JSON.stringify(controls.target.toArray());}
 connector.replaceChildren();const connected=basis==='model'&&layer==='species'&&scale!=='mouth'&&scale!=='biofilm';connector.toggleAttribute('hidden',!connected);
 if(connected){
  const p=new T.Vector3(4.55,(PRESETS[health].margin+PRESETS[health].attachment)/2,0).project(camera),r=lens.getBoundingClientRect(),s=stage.getBoundingClientRect(),x=fallback?stage.clientWidth*.64:(p.x+1)/2*stage.clientWidth,y=fallback?stage.clientHeight*.45:(-p.y+1)/2*canvas.clientHeight+62;
  if(Number.isFinite(x)&&Number.isFinite(y)){  connector.setAttribute('viewBox',`0 0 ${stage.clientWidth} ${stage.clientHeight}`);const line=document.createElementNS('http://www.w3.org/2000/svg','path');line.setAttribute('d',`M${x},${y} L${r.left-s.left+9},${r.top-s.top+9}`);line.setAttribute('stroke','#a4bda8');line.setAttribute('stroke-opacity','.4');line.setAttribute('stroke-dasharray','3 5');line.setAttribute('fill','none');connector.append(line);}
 }
 leaders.replaceChildren();
 labels.hidden=!showLabels||scale==='biofilm';
 if(scale==='mouth'){
  for(const [id,b]of labelElements){b.hidden=id!=='plaque';if(id==='plaque'){b.textContent='Selected FDI 36 →';const p=new T.Vector3(24.4,7,9.2).project(camera);b.style.left=((p.x+1)/2*stage.clientWidth)+'px';b.style.top=((-p.y+1)/2*canvas.clientHeight+62+24)+'px';}}
  return;
 }
 const h=canvas.clientHeight||stage.clientHeight;
 const bounds=modelBox(),projected:Array<T.Vector3>=[];
 for(const x of [bounds.min.x,bounds.max.x])for(const y of [bounds.min.y,bounds.max.y])for(const z of [bounds.min.z,bounds.max.z])projected.push(new T.Vector3(x,y,z).project(camera));
 const xs=projected.map(p=>(p.x+1)/2*stage.clientWidth).filter(Number.isFinite),minX=xs.length?Math.min(...xs):0,maxX=xs.length?Math.max(...xs):stage.clientWidth;
 stage.dataset.modelBounds=JSON.stringify({left:minX,right:maxX});
 leaders.setAttribute('viewBox',`0 0 ${stage.clientWidth} ${stage.clientHeight}`);
 const left:Array<{b:HTMLButtonElement;x:number;y:number;ax:number;ay:number}>=[],right:Array<{b:HTMLButtonElement;x:number;y:number;ax:number;ay:number}>=[];
 for(const [id,b]of labelElements){b.hidden=(!cutaway&&['dentin','pulp','cementum','plaque','lumen','epithelium','pdl'].includes(id))||(['cheek','tongue'].includes(id)&&scale!=='tooth');if(b.hidden)continue;b.textContent=({pdl:'Ligament',epithelium:'Junctional lining',supragingival:'Supragingival plaque',cheek:'Cheek side · buccal',tongue:'Tongue side · lingual'} as Partial<Record<PocketId,string>>)[id]||names[id];b.dataset.selected=String(id===selected);
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
 if(scale==='mouth'){inspect('FDI 36 in the mouth','The mandibular left first permanent molar.','The marker locates the selected tooth within authored upper and lower dental arches, palate and tongue context.','This is an illustrative placement, not a measured registration. Oral tissues and tooth morphology await expert anatomical review.','Authored oral context · no physical scale','anatomy');return;}
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
function selectTissue(id:PocketId){
 if(scale==='mouth'&&id==='plaque'){setScale('tooth');return;}
 selected=id;get<HTMLSelectElement>('tissue').value=id;showInspector();if(id==='plaque'&&layer==='species'){inspectSpecies();if(scale!=='biofilm')setScale('biofilm');}
 else if(id==='supragingival')inspect(names[id],...descriptions[id],'This viewer assigns no species positions to supragingival plaque. The selected spatial observations concern attached subgingival biofilm in advanced periodontitis and are not transferred here.','Authored compartment · species locations unresolved','plaque');
 else inspect(names[id],...descriptions[id],SOURCES.anatomy.method+' '+SOURCES.anatomy.limit,'Authored anatomy · expert review pending','anatomy');
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
function renderLens(){
 const visible=layer==='species'&&scale!=='mouth';lens.hidden=!visible;stage.classList.toggle('detail',scale==='biofilm');stage.classList.toggle('with-lens',visible&&scale!=='biofilm');
 get('detail').textContent=scale==='biofilm'?'← Return to pocket':'Explore biofilm detail →';
 const run=saved.cases.find(c=>c.id===health)!;
 if(basis==='evidence'){
  get('lens-title').textContent='Across attached plaque';get('lens-note').textContent=health==='healthy'?'Healthy-site locations unresolved; disease-specimen observations are not transferred.':'Qualitative human observations · schematic glyphs, not counts or measured positions.';
  get('lens-content').innerHTML=evidenceSVG();get('lens-axis').textContent='Across-biofilm organization ≠ bands down the pocket. No physical scale.';
 }else{
  get('lens-title').textContent='Assumed model coordinates';get('lens-note').textContent='Original 32 × 72 grid · fixed saved endpoint · no anatomical projection.';
  const c=document.createElement('canvas');c.width=run.nx;c.height=run.ny;c.setAttribute('aria-label','Saved '+health+' '+taxon+' density grid, original coordinates');c.setAttribute('role','img');c.dataset.case=run.id;c.dataset.taxon=taxon;c.dataset.sourceSha256=saved.source.sha256;c.dataset.frame=String(run.frame_index);
  const ctx=c.getContext('2d')!,pixels=ctx.createImageData(run.nx,run.ny);pixels.data.set(gridRGBA(run,taxon));ctx.putImageData(pixels,0,0);
  const image=document.createElement('img');image.id='grid-image';image.src=c.toDataURL('image/png');image.alt=c.getAttribute('aria-label')!;image.width=run.nx;image.height=run.ny;for(const [key,value]of Object.entries(c.dataset))image.dataset[key]=value;get('lens-content').replaceChildren(image);
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
 stage.dataset.health=health;stage.dataset.scale=scale;stage.dataset.layer=layer;stage.dataset.basis=basis;stage.dataset.taxon=taxon;stage.dataset.cutaway=String(cutaway);stage.dataset.host='human';stage.dataset.anatomy='Blender-authored generic model; expert review pending';
 for(const [attr,value]of [['health',health],['scale',scale==='biofilm'?'pocket':scale],['layer',layer],['basis',basis]])for(const b of document.querySelectorAll<HTMLButtonElement>(`[data-${attr}]`))b.setAttribute('aria-pressed',String(b.dataset[attr]===value));
 get('labels').setAttribute('aria-pressed',String(showLabels));get('cutaway').setAttribute('aria-pressed',String(cutaway));get<HTMLButtonElement>('cutaway').disabled=scale==='mouth'||scale==='biofilm'||fallback;
 document.querySelector<HTMLButtonElement>('[data-scale="mouth"]')!.disabled=fallback;
 get('scene-title').textContent=scale==='mouth'?'Oral cavity · selected FDI 36':scale==='biofilm'?'Biofilm detail · '+PRESETS[health].title:PRESETS[health].title+(cutaway?' · quarter cutaway':' · assembled');
 get('scene-subtitle').textContent=scale==='mouth'?'Authored arches, palate and tongue · no measured tooth registration':'Illustrative landmarks · no physical scale';
 get('orientation').hidden=scale==='biofilm';get('structure-controls').hidden=layer!=='structure';get('species-controls').hidden=layer!=='species';document.querySelector<HTMLElement>('.basis-controls')!.hidden=layer!=='species';get('environment-controls').hidden=layer!=='environment';get('environment-content').hidden=layer!=='environment';
 get<HTMLSelectElement>('taxon').value=taxon;
 for(const b of get('legend').querySelectorAll<HTMLButtonElement>('button'))b.setAttribute('aria-pressed',String(b.dataset.taxon===taxon));
 const note=get('scene-note');
 if(layer==='species'){
  note.textContent=basis==='model'?'Assumed model coordinates · Dominant color identifies the largest modeled density in a bin; other taxa may overlap. Saved fields have unresolved numerical and biological limitations.':health==='healthy'?'Location unresolved for healthy anatomy. The selected spatial evidence concerns advanced periodontitis; no disease-specimen positions are copied.':'Source-scoped observations across attached biofilm. These schematic arrangements do not establish universal pocket-depth positions, cell contacts or abundance.';
  inspectSpecies();
 }else if(layer==='environment'){note.textContent='Source observations remain separate from site-specific values. No continuous oxygen/pH field is invented.';renderEnvironment();}
 else {note.textContent='Margin, epithelial attachment and bone crest are separate authored presets. Geometry remains uncalibrated; expert anatomical review is pending.';inspect(names[selected],...descriptions[selected],SOURCES.anatomy.method+' '+SOURCES.anatomy.limit,'Authored anatomy · expert review pending','anatomy');}
 renderLens();resize();highlight();
}
function setScale(next:Scale){
 remember();const previous=scale;scale=next;
 if(next!=='biofilm'){geometry();view();}else if(previous==='mouth'){geometry();view();}
 refresh();resize();
}
function resize(){if(renderer){const h=stage.classList.contains('with-lens')?stage.clientHeight-247:stage.clientHeight-62;renderer.setPixelRatio(Math.min(devicePixelRatio,1.5));renderer.setSize(stage.clientWidth,h,false);camera.aspect=stage.clientWidth/h;camera.updateProjectionMatrix();composer?.setPixelRatio(renderer.getPixelRatio());composer?.setSize(stage.clientWidth,h);}updateLabels();dirty=true;}
for(const b of document.querySelectorAll<HTMLButtonElement>('[data-health]'))b.onclick=()=>{health=b.dataset.health as Health;geometry();refresh();};
for(const b of document.querySelectorAll<HTMLButtonElement>('[data-scale]'))b.onclick=()=>setScale(b.dataset.scale as Scale);
for(const b of document.querySelectorAll<HTMLButtonElement>('[data-layer]'))b.onclick=()=>{layer=b.dataset.layer as Layer;if(layer!=='species'&&scale==='biofilm')setScale('pocket');refresh();};
for(const b of document.querySelectorAll<HTMLButtonElement>('[data-basis]'))b.onclick=()=>{basis=b.dataset.basis as Basis;refresh();};
const select=get<HTMLSelectElement>('taxon');
const tissueSelect=get<HTMLSelectElement>('tissue');for(const id of ['enamel','dentin','pulp','cementum','pdl','gingiva','bone','supragingival','plaque','lumen','epithelium','connective','vessels','nerve','cheek','tongue'] as PocketId[]){const option=document.createElement('option');option.value=id;option.textContent=names[id];tissueSelect.append(option);}tissueSelect.value='plaque';tissueSelect.onchange=()=>selectTissue(tissueSelect.value as PocketId);
for(const t of TAXA){const option=document.createElement('option');option.value=t.id;option.textContent=t.name;select.append(option);const b=document.createElement('button'),swatch=document.createElement('i');swatch.style.background=t.color;b.dataset.taxon=t.id;b.setAttribute('aria-pressed','false');b.append(swatch,document.createTextNode(t.short));b.onclick=()=>{taxon=t.id;showInspector();refresh();};get('legend').append(b);}
select.onchange=()=>{taxon=select.value;showInspector();refresh();};
get('detail').onclick=()=>setScale(scale==='biofilm'?'pocket':'biofilm');
get('labels').onclick=()=>{showLabels=!showLabels;refresh();};get('cutaway').onclick=()=>{cutaway=!cutaway;geometry();refresh();};get('reset').onclick=()=>{if(scale==='biofilm')setScale('pocket');remembered.delete(scale);view(true);};
get('close-inspector').onclick=()=>{inspector=false;get('inspector').hidden=true;get('workspace').classList.add('inspector-closed');get('show-inspector').hidden=false;resize();};get('show-inspector').onclick=showInspector;
get<HTMLSelectElement>('compartment').onchange=renderEnvironment;
get('provenance').textContent='Source: '+saved.source.file+' · SHA-256 '+saved.source.sha256+'. '+saved.source.precision+' Both cases show their fixed endpoint; model time is not physical elapsed time.';
get('render-check').onclick=()=>{if(!renderer){get('performance').textContent='2D fallback active; WebGL performance not measured.';return;}if(scale==='biofilm')setScale('pocket');measuring=true;measureStart=performance.now();frameCount=0;renderMS=[];get('performance').textContent='Measuring 60 rendered frames…';dirty=true;};
canvas.addEventListener('keydown',e=>{if(!controls)return;if(e.key==='Home'){e.preventDefault();view(true);return;}if(!e.key.startsWith('Arrow'))return;e.preventDefault();const offset=camera.position.clone().sub(controls.target),s=new T.Spherical().setFromVector3(offset);s.theta+=e.key==='ArrowLeft'?-.10:e.key==='ArrowRight'?.10:0;s.phi=T.MathUtils.clamp(s.phi+(e.key==='ArrowUp'?-.08:e.key==='ArrowDown'?.08:0),.1,Math.PI-.1);camera.position.copy(controls.target).add(new T.Vector3().setFromSpherical(s));controls.update();dirty=true;updateLabels();});
let down=[0,0];canvas.addEventListener('pointerdown',e=>{down=[e.clientX,e.clientY];});canvas.addEventListener('pointerup',e=>{if(!renderer||Math.hypot(e.clientX-down[0],e.clientY-down[1])>5)return;const r=canvas.getBoundingClientRect(),pointer=new T.Vector2((e.clientX-r.left)/r.width*2-1,-(e.clientY-r.top)/r.height*2+1),ray=new T.Raycaster();ray.setFromCamera(pointer,camera);const hits=ray.intersectObjects(display).filter(h=>h.object.visible),hit=hits.find(h=>!h.object.userData.context)||hits[0];if(hit)selectTissue(hit.object.userData.tissue);});
function fallbackScene(reason:string){fallback=true;stage.dataset.renderer='2d-fallback';canvas.hidden=true;get('fallback').hidden=false;status(reason+' The 2D anatomy and all evidence/model controls remain available.');updateFallback();}
try{
 if(new URLSearchParams(location.search).get('render')==='2d')throw Error('2D fallback requested');
 renderer=new T.WebGLRenderer({canvas,antialias:true,alpha:true,powerPreference:'high-performance'});renderer.setClearColor('#0c171a',0);renderer.toneMapping=T.NeutralToneMapping;renderer.toneMappingExposure=.74;renderer.shadowMap.enabled=true;renderer.shadowMap.type=T.PCFShadowMap;
 mats=materials();boneCapMap=boneSectionTexture();fiberCapMap=connectiveSectionTexture();for(const m of Object.values(mats))m.envMapIntensity=.28;mats.enamel.envMapIntensity=.42;scene.environment=environment(renderer);
 // Restrained three-point lighting: warm soft key, cool low fill, neutral rim for silhouette separation.
 const key=new T.DirectionalLight('#ffe9d6',2.1);key.position.set(-15,24,18);key.castShadow=true;key.shadow.mapSize.set(2048,2048);key.shadow.camera.left=-25;key.shadow.camera.right=25;key.shadow.camera.top=25;key.shadow.camera.bottom=-25;key.shadow.bias=-.0002;key.shadow.normalBias=.02;key.shadow.radius=3;
 const fill=new T.DirectionalLight('#b4cfe0',.55);fill.position.set(18,2,16);const rim=new T.DirectionalLight('#e8eef2',1.25);rim.position.set(-8,12,-18);scene.add(key,fill,rim,new T.HemisphereLight('#e4e2dc','#1a2124',.22));
 controls=new OrbitControls(camera,canvas);controls.enableDamping=false;controls.enablePan=false;controls.minDistance=20;controls.maxDistance=145;controls.addEventListener('change',()=>{dirty=true;updateLabels();});
 composer=new EffectComposer(renderer);composer.addPass(new RenderPass(scene,camera));
 ao=new SSAOPass(scene,camera,512,512,16);ao.kernelRadius=.8;ao.minDistance=.0002;ao.maxDistance=.025;composer.addPass(ao);composer.addPass(new OutputPass());
 stage.dataset.renderEffects='Contact shadows · physical materials · filmic tone mapping';
 stage.dataset.renderer='webgl2';geometry();resize();view(true);status('Local 3D scene ready · source-scoped evidence · expert anatomical review pending.');
 canvas.addEventListener('webglcontextlost',e=>{e.preventDefault();fallbackScene('WebGL context was lost.');refresh();});
}catch(error){fallbackScene((error as Error).message+'.');}
refresh();resize();if(controls)view(true);document.body.dataset.ready='true';
new ResizeObserver(resize).observe(stage);
function tick(){requestAnimationFrame(tick);if(renderer&&!fallback&&(dirty||measuring)&&scale!=='biofilm'){
 const start=performance.now();composer?composer.render():renderer.render(scene,camera);dirty=false;updateLabels();
 if(measuring){renderMS.push(performance.now()-start);frameCount++;if(frameCount>=60){measuring=false;const wall=performance.now()-measureStart,sorted=renderMS.slice().sort((a,b)=>a-b),fps=frameCount*1000/wall;stage.dataset.measuredFps=fps.toFixed(1);stage.dataset.renderMedianMs=sorted[Math.floor(sorted.length/2)].toFixed(2);get('performance').textContent=` ${fps.toFixed(1)} frames/s over ${frameCount} frames · median CPU submission ${sorted[Math.floor(sorted.length/2)].toFixed(2)} ms · ${Math.round(stage.clientWidth)} × ${Math.round(stage.clientHeight)} CSS pixels. This measures display execution, not biological accuracy.`;}}
 }}tick();
