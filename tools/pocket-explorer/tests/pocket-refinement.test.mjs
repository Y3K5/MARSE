import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
globalThis.window={};
const {buildPocket,SEGMENT,ANCHORS}=await import('../dist/pocket-core/pocket-geometry.mjs');
const {endSection}=await import('../dist/pocket-core/pocket-section.mjs');
const {PocketViewCache}=await import('../dist/pocket-core/pocket-view-cache.mjs');
const {visibleLabels,sceneLabelIds,sceneCopy,additionalLabelTargets}=await import('../dist/pocket-core/pocket-layout.mjs');
const asset=JSON.parse(readFileSync(new URL('../assets/pocket-organic/anatomy.json',import.meta.url)));

test('every tissue closes at both specimen-end planes',()=>{
 for(const health of ['healthy','periodontitis'])for(const part of buildPocket(health)){
  const section=endSection(part.geometry,...SEGMENT.ends);
  assert.equal(section.openChains,0,health+' '+part.name);
  const p=section.surface.getAttribute('position');
  for(let i=0;i<p.count;i++)assert.ok(p.getX(i)>=SEGMENT.ends[0]-1e-4&&p.getX(i)<=SEGMENT.ends[1]+1e-4);
  section.surface.dispose();section.cap.dispose();part.geometry.dispose();
 }
});
test('packed format v2 decodes 32-bit indices for every part',()=>{
 const packed=JSON.parse(readFileSync(new URL('../src/pocket-organic-input.json',import.meta.url)));
 assert.equal(packed.format,'marse.pocket-organic-packed/2');assert.match(packed.precision,/Uint32/);
 for(const health of ['healthy','periodontitis'])for(const part of buildPocket(health)){
  assert.ok(part.geometry.index.array instanceof Uint32Array);part.geometry.dispose();
 }
});
test('both presets stay within the 600k Pocket-view triangle budget and reuse prepared sections',()=>{
 const cache=new PocketViewCache();
 try {for(const health of ['healthy','periodontitis'])for(const cut of [false,true]){
  const view=cache.get(health,cut),built=cache.builds;
  assert.equal(view.openChains,0,health+' '+cut);
  assert.ok(view.triangles<=600000,health+' cutaway='+cut+' uses '+view.triangles+' triangles');
  assert.equal(cache.get(health,cut),view);assert.equal(cache.builds,built);
  assert.equal(cache.get(health,cut).parts[0].geometry,view.parts[0].geometry);
 }}finally{cache.dispose();}
});
test('small-stage labels remain compact while a selected tissue stays reachable',()=>{
 for(const id of sceneLabelIds){assert.ok(ANCHORS[id]??additionalLabelTargets[id],id+' needs a label target');const labels=visibleLabels(390,true,'tooth',id);assert.ok(labels.length<=4,id);assert.ok(labels.includes(id),id);}
 assert.ok(visibleLabels(768,true,'pocket','plaque').length>4);
 assert.ok(!visibleLabels(390,false,'pocket','plaque').includes('lumen'));
 assert.match(sceneCopy('tooth','healthy',false).title,/FDI 35–37/);
 assert.match(sceneCopy('pocket','periodontitis',true).title,/pocket/);
});
test('neighbor crowns differ and selected-site landmarks stay fixed',()=>{
 const f=asset.authored_features;
 assert.equal(f.context_crowns['FDI 35'].cusps,3);
 assert.equal(f.context_crowns['FDI 37'].cusps,4);assert.equal(f.context_crowns['FDI 37'].grooves,'cross');
 assert.deepEqual(f.pocket_landmarks,{healthy:{margin:.35,attachment:-1.65},periodontitis:{margin:-.15,attachment:-6.2}});
 const frozen=JSON.parse(readFileSync(new URL('../fixtures/authored-landmarks.json',import.meta.url)));
 assert.deepEqual(f.distal_crest_samples.healthy,frozen.healthy_crest_samples);
 assert.deepEqual(f.distal_crest_samples.periodontitis.filter(s=>s.x<=4.6),frozen.periodontitis_selected_crest_samples);
});
