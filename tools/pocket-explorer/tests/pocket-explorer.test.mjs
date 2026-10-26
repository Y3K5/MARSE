import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {TAXA,PRESETS,observation,gridRGBA,validateSaved} from '../dist/pocket-core/pocket-science.mjs';
import {sectionGeometry} from '../dist/core.mjs';
globalThis.window={};
const {buildPocket}=await import('../dist/pocket-core/pocket-geometry.mjs');
const raw=readFileSync(new URL('../src/pocket-input.json',import.meta.url)),input=JSON.parse(raw),baseline=JSON.parse(readFileSync(new URL('../fixtures/endpoints.json',import.meta.url)));
test('preserved endpoints retain their reviewed identities without importing full trajectories',()=>{
 assert.equal(createHash('sha256').update(raw).digest('hex'),baseline.input_sha256);
 assert.deepEqual(input.source,baseline.source);assert.equal(validateSaved(input),input);
 for(const c of input.cases){const original=baseline.cases.find(s=>s.id===c.id);for(const key of ['frame_index','model_time','nx','ny'])assert.equal(c[key],original[key]);assert.equal(createHash('sha256').update(JSON.stringify(c.biomass)).digest('hex'),original.biomass_sha256);}
});
test('every taxon can be rendered without changing the input, and dominant identity preserves overlap',()=>{
 const before=JSON.stringify(input);
 for(const c of input.cases){const signatures=new Set();for(const t of TAXA){const pixels=gridRGBA(c,t.id);assert.equal(pixels.length,32*72*4);signatures.add(createHash('sha256').update(pixels).digest('hex'));}assert.equal(signatures.size,9);}
 const c={...input.cases[0],biomass:TAXA.map(()=>Array(2304).fill(0))};c.biomass[0][0]=2;c.biomass[1][0]=1;
 assert.deepEqual(Array.from(gridRGBA(c,'all').slice(0,3)),[126,207,192]);assert.equal(c.biomass[1][0],1);assert.equal(JSON.stringify(input),before);
});
test('unknown locations and healthy host-state evidence stay unresolved',()=>{
 for(const t of TAXA)assert.equal(observation(t.id,'healthy').resolved,false,t.id);
 for(const id of ['ss','vp','td'])assert.equal(observation(id,'periodontitis').resolved,false,id);
 assert.equal(TAXA.find(t=>t.id==='ac').level,'genus/group');assert.equal(observation('fn','periodontitis').region,'intermediate');
});
test('foreign hosts, reordered identities and malformed grids are rejected',()=>{
 for(const mutate of [d=>d.host='canine',d=>d.source.sha256='unknown',d=>d.taxa.reverse(),d=>d.cases.pop(),d=>d.cases[0].biomass[0][0]=-1,d=>d.cases[0].biomass[0].pop(),d=>d.units='mol/m3',d=>d.coordinates.order='x*ny+y',d=>d.coordinates.physical_registration='measured',d=>d.source.precision='']){const copy=structuredClone(input);mutate(copy);assert.throws(()=>validateSaved(copy));}
});
test('both presets have separate anatomical landmarks and finite closed sampled cutaways',()=>{
 assert.notEqual(PRESETS.healthy.margin,PRESETS.periodontitis.margin);assert.notEqual(PRESETS.healthy.attachment,PRESETS.periodontitis.attachment);assert.notEqual(PRESETS.healthy.crest,PRESETS.periodontitis.crest);
 for(const health of ['healthy','periodontitis']){
  const parts=buildPocket(health),ids=new Set(parts.map(p=>p.id));for(const id of ['enamel','dentin','pulp','cementum','plaque','lumen','epithelium','connective','pdl','bone'])assert.ok(ids.has(id));
  for(const p of parts){assert.ok(Array.from(p.geometry.getAttribute('position').array).every(Number.isFinite),p.name);for(const cut of [-1.1,0,1.3]){const section=sectionGeometry(p.geometry,cut);assert.equal(section.openChains,0,health+' '+p.name+' '+cut);section.surface.dispose();section.cap.dispose();}p.geometry.dispose();}
 }
});

test('Blender tissues are closed and the dentin crown and roots form one connected body',()=>{
 const asset=JSON.parse(readFileSync(new URL('../assets/pocket-organic/anatomy.json',import.meta.url)));
 assert.equal(asset.calibration,null);assert.equal(asset.review,'pending');
 for(const rec of asset.topology){assert.equal(rec.nonmanifold_edges,0,rec.case+' '+rec.part);if(['Continuous crown and root dentin','Continuous pulp chamber horns and tapered canals'].includes(rec.part))assert.equal(rec.components,1,rec.case+' '+rec.part);}
});

test('packed browser meshes preserve the canonical Blender vertices and indices',()=>{
 const source=JSON.parse(readFileSync(new URL('../assets/pocket-organic/anatomy.json',import.meta.url)));
 for(const health of ['healthy','periodontitis']){const actual=buildPocket(health);source.cases[health].forEach((part,i)=>{assert.equal(actual[i].id,part.id);assert.deepEqual(actual[i].geometry.getAttribute('position').array,new Float32Array(part.positions));assert.deepEqual(actual[i].geometry.index.array,new Uint32Array(part.indices));actual[i].geometry.dispose();});}
});

const {quarterSection}=await import('../dist/pocket-core/pocket-section.mjs');
const THREE=await import('three');
test('quarter cutaway keeps exactly three quarters of a cube and closes both exposed planes',()=>{
 const original=new THREE.BoxGeometry(2,2,2),s=quarterSection(original);assert.equal(s.openChains,0);
 const volume=g=>{const p=g.getAttribute('position');let sum=0;for(let i=0;i<p.count;i+=3){const a=new THREE.Vector3().fromBufferAttribute(p,i),b=new THREE.Vector3().fromBufferAttribute(p,i+1),c=new THREE.Vector3().fromBufferAttribute(p,i+2);sum+=a.dot(b.cross(c))/6;}return sum;};
 assert.ok(Math.abs(volume(s.surface)+volume(s.cap)-6)<1e-10);
 const cp=s.cap.getAttribute('position');for(let i=0;i<cp.count;i++)assert.ok(Math.abs(cp.getX(i))<1e-6||Math.abs(cp.getZ(i))<1e-6);
 original.dispose();s.surface.dispose();s.cap.dispose();
});
test('every organic tissue has closed quarter-cut construction and finite normals',()=>{
 for(const health of ['healthy','periodontitis'])for(const p of buildPocket(health)){
  const s=quarterSection(p.geometry);assert.equal(s.openChains,0,health+' '+p.name);
  for(const g of [s.surface,s.cap]){assert.ok(Array.from(g.getAttribute('position').array).every(Number.isFinite));assert.ok(Array.from(g.getAttribute('normal').array).every(Number.isFinite));g.dispose();}p.geometry.dispose();
 }
});
