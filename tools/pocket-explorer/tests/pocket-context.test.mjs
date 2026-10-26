import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import * as T from 'three';
globalThis.window={};
const {oralContext,facialContext,SEGMENT_TO_MOUTH,FDI36_POSE,neuralContext}=await import('../dist/pocket-core/pocket-context.mjs');
const {CROWN_REFERENCE,NERVE_STEPS,PARAMETERS}=await import('../dist/pocket-core/pocket-reference.mjs');
const {buildPocket}=await import('../dist/pocket-core/pocket-geometry.mjs');
const {SEQUENCES,stageAt,sequenceSVG}=await import('../dist/pocket-core/pocket-sequences.mjs');
const volume=g=>{const p=g.getAttribute('position'),ids=g.index;let v=0;for(let i=0;i<ids.count;i+=3){const a=new T.Vector3().fromBufferAttribute(p,ids.getX(i)),b=new T.Vector3().fromBufferAttribute(p,ids.getX(i+1)),c=new T.Vector3().fromBufferAttribute(p,ids.getX(i+2));v+=a.dot(b.cross(c))/6;}return v;};
test('context uses 27 named crowns plus the preserved FDI36; all context crowns face outward',()=>{
 const parts=oralContext(),crowns=parts.filter(p=>p.name.startsWith('FDI '));
 assert.equal(crowns.length,27);assert.equal(new Set(crowns.map(p=>p.name.split(' · ')[0])).size,27);
 assert.ok(!crowns.some(p=>p.name.startsWith('FDI 36')));
 for(const p of parts){assert.ok(Array.from(p.geometry.getAttribute('position').array).every(Number.isFinite));if(crowns.includes(p))assert.ok(volume(p.geometry)>0,p.name);p.geometry.dispose();}
 const origin=new T.Vector3().applyMatrix4(SEGMENT_TO_MOUTH);assert.ok(origin.distanceTo(new T.Vector3(...FDI36_POSE.position))<1e-12);
 const buccal=new T.Vector3(0,0,1).transformDirection(SEGMENT_TO_MOUTH);assert.ok(buccal.x>0,'buccal faces subject left');
});
test('chosen crown reference affects cloned presentation geometry, preserving source meshes',()=>{
 const parts=buildPocket('periodontitis'),p=parts.find(p=>p.id==='enamel'&&!p.name.includes('Neighbouring'));
 const before=p.geometry.getAttribute('position').array.slice(),g=p.geometry.clone().applyMatrix4(new T.Matrix4().makeScale(...CROWN_REFERENCE.scale));g.computeBoundingBox();const size=g.boundingBox.getSize(new T.Vector3());
 assert.ok(Math.abs(size.x-CROWN_REFERENCE.md)<2e-6);assert.ok(Math.abs(size.z-CROWN_REFERENCE.bl)<2e-6);assert.deepEqual(p.geometry.getAttribute('position').array,before);
 assert.ok(PARAMETERS.some(p=>p.classification==='Unresolved'));g.dispose();for(const part of parts)part.geometry.dispose();
});
test('face is a pinned CC0 graphical asset with explicit authored registration',()=>{
 const face=JSON.parse(readFileSync(new URL('../src/pocket-face-input.json',import.meta.url)));
 assert.equal(face.license,'CC0-1.0');assert.equal(face.commit,'a8bc2d54ff0ac92e78ff71431b1023eda42bf482');assert.match(face.registration.status,/not anatomically measured/);
 assert.equal(face.triangles,face.indices.length/3);assert.ok(face.positions.every(Number.isFinite));assert.ok(face.indices.every(i=>Number.isInteger(i)&&i>=0&&i*3<face.positions.length));
 for(const part of facialContext()){assert.ok(part.geometry.getAttribute('normal'));part.geometry.dispose();}
});
test('storyboards replay deterministically without changing presets or saved fields',()=>{
 const saved=readFileSync(new URL('../src/pocket-input.json',import.meta.url)),before=JSON.stringify(SEQUENCES);
 for(const seq of ['assembly','response'])for(let t=0;t<=5;t+=.25){const first=sequenceSVG(seq,t);assert.equal(first,sequenceSVG(seq,t));assert.match(first,/anonymous illustrative storyboard/);assert.ok(stageAt(seq,t).source);}
 assert.equal(stageAt('assembly',-2).index,0);assert.equal(stageAt('assembly',9).index,5);assert.throws(()=>stageAt('assembly',NaN));assert.equal(JSON.stringify(SEQUENCES),before);assert.deepEqual(readFileSync(new URL('../src/pocket-input.json',import.meta.url)),saved);
});
test('sensory connections distinguish dental from mental supply and vessels remain separate',()=>{
 assert.equal(NERVE_STEPS.length,6);assert.match(NERVE_STEPS.map(s=>s.text).join(' '),/mental/i);
 const parts=neuralContext(true);for(const style of ['nerve','artery','vein'])assert.ok(parts.some(p=>p.appearance===style));
 assert.ok(parts.some(p=>p.name.startsWith('Mental branch')));for(const part of parts)part.geometry.dispose();
});
