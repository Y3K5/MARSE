// Local geometry/software controls. These do not establish biological accuracy.
import {readFile,writeFile,mkdtemp} from 'node:fs/promises';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {pathToFileURL} from 'node:url';
import assert from 'node:assert/strict';
const home=process.argv[2]||new URL('.',import.meta.url).pathname;
const stage=process.env.MARSE_TEST_STAGE||home;
const threeURL=pathToFileURL(join(home,'vendor/three.module.js')).href;
const T=await import(threeURL),tmp=await mkdtemp(join(tmpdir(),'marse-geometry-'));
await writeFile(join(tmp,'flow.mjs'),(await readFile(join(stage,'flow-dive.js'),'utf8')).replace("from 'three'",`from '${threeURL}'`));
const {FlowDive}=await import(pathToFileURL(join(tmp,'flow.mjs')));
globalThis.matchMedia=()=>({matches:true});globalThis.document={documentElement:{dataset:{}}};
const checks=[];
function setup(geometry,record=null){const group=new T.Group(),surface=new T.Mesh(geometry,new T.MeshStandardMaterial());surface.userData.cellBoundary=true;group.add(surface);const camera=new T.PerspectiveCamera(38,1.7,.01,100);camera.position.set(0,0,12);const controls={target:new T.Vector3(),enabled:true,update(){camera.lookAt(this.target);}};return new FlowDive({group,camera,controls,record,organism:'software-fixture'});}
const pts=Array.from({length:65},(_,i)=>new T.Vector3((i/64-.5)*8,Math.sin(i/64*Math.PI*5)*.5,Math.cos(i/64*Math.PI*5)*.5));
for(const [name,geometry]of [['thin spirochete',new T.TubeGeometry(new T.CatmullRomCurve3(pts),160,.14,24)],['coccus',new T.SphereGeometry(.4,16,12)],['rod',new T.CapsuleGeometry(.4,1.3,6,12)]]){
 const f=setup(geometry);assert(f.patchSamples.length>20);assert(f.patchRadius<=f.localRadius);assert(f.info().surfaceResidual<1e-10);for(const sample of f.patchSamples)assert(sample.n.length()>.99);const ids=[f.group.uuid,f.root.uuid,f.protein.uuid];for(const level of [0,1,2,0]){f.setLevel(level);f.update();assert.equal(f.level,level);assert.deepEqual([f.group.uuid,f.root.uuid,f.protein.uuid],ids);assert(f.camera.near>0);}checks.push(name+': connected samples, normal offsets, residual and persistent identity');
}
// Two close parallel surfaces in one mesh, with disconnected topology.
const g=new T.BufferGeometry(),positions=[];for(const z of [0,.001])positions.push(-1,-1,z,1,-1,z,1,1,z,-1,-1,z,1,1,z,-1,1,z);g.setAttribute('position',new T.Float32BufferAttribute(positions,3));g.computeVertexNormals();const disconnected=setup(g);const selected=disconnected.patchSamples.map(s=>s.triangle.every(i=>disconnected.graph.vertices[i].p.z>.0005));assert(selected.every(Boolean)||selected.every(v=>!v));checks.push('disconnected near surfaces cannot enter one region');
const record={id:'constructed-software-control',sha256:'fixture',chains:{A:[{p:[0,0,0],r:'1'},{p:[3,0,0],r:'2'},{p:[3,4,0],r:'4'},{p:[3,7,0],r:'5'}]}};
const protein=setup(new T.SphereGeometry(1,16,12),record);assert.equal(protein.sourcePoints.length,4);assert.equal(protein.backbone.children[0].children.length,2);assert(protein.info().coordinatePreservationMaxPairError<1e-12);assert(protein.info().renderedCoordinateError<1e-4);assert.equal(protein.protein.scale.x,protein.protein.scale.y);assert.equal(protein.protein.scale.y,protein.protein.scale.z);checks.push('rigid source frame, uniform scaling, instance coordinates and sequence-gap breaks');
const sourceCode=await readFile(join(stage,'source-protein.js'),'utf8');await writeFile(join(tmp,'source.mjs'),sourceCode);const {loadSourceProtein}=await import(pathToFileURL(join(tmp,'source.mjs')));
globalThis.fetch=async()=>({ok:false});await assert.rejects(loadSourceProtein('missing-meta'),/metadata unavailable/);
globalThis.fetch=async url=>url.endsWith('.json')?{ok:true,json:async()=>({format:'pdb'})}:{ok:false};await assert.rejects(loadSourceProtein('missing-file'),/coordinate file unavailable/);
globalThis.fetch=async url=>url.endsWith('.json')?{ok:true,json:async()=>({sha256:'0'.repeat(64),format:'pdb'})}:{ok:true,arrayBuffer:async()=>new TextEncoder().encode('constructed control').buffer};await assert.rejects(loadSourceProtein('wrong-hash'),/checksum mismatch/);checks.push('unavailable metadata, unavailable coordinates and checksum errors refuse before parsing');
console.log(JSON.stringify({checks,passed:checks.length,scope:'software and geometry only; no biological validation'},null,2));
