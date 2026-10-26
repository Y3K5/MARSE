import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import * as T from 'three';
import {atlasContext,ATLAS_BASIS,ATLAS_OPEN,ATLAS_FDI36} from '../dist/pocket-core/pocket-atlas.mjs';
const input=JSON.parse(await readFile(new URL('../src/pocket-atlas-input.json',import.meta.url)));

test('atlas keeps 28 named teeth, exact source identities and separate source/pose claims',()=>{
 const teeth=input.parts.filter(p=>p.fdi);
 assert.equal(teeth.length,28);assert.equal(new Set(teeth.map(p=>p.fdi)).size,28);
 assert.deepEqual(teeth.find(p=>p.fdi===36).concept,'FMA55704');
 assert.equal(teeth.find(p=>p.fdi===36).element,'FJ1254');
 for(const p of input.parts)assert.match(p.source_sha256,/^[a-f0-9]{64}$/);
 assert.equal(input.license,'CC-BY-4.0');assert.match(input.specimen,/not an average/);
});
test('source positions survive the common basis change and authored jaw pose',()=>{
 const original=JSON.stringify(input);const parts=atlasContext('face');
 assert.equal(ATLAS_BASIS.determinant(),1);assert.ok(Math.abs(ATLAS_OPEN.determinant()-1)<1e-14);
 for(const p of parts){
  const source=input.parts.find(q=>q.element===p.atlasIdentity.element);
  const buffer=Buffer.from(source.positions,'base64');
  const a=p.geometry.getAttribute('position');
  for(let i=0;i<a.count;i++){
   const expected=new T.Vector3(buffer.readFloatLE(i*12),buffer.readFloatLE(i*12+4),buffer.readFloatLE(i*12+8)).applyMatrix4(ATLAS_BASIS);
   if(source.jaw==='lower')expected.applyMatrix4(ATLAS_OPEN);
   assert.ok(expected.distanceTo(new T.Vector3().fromBufferAttribute(a,i))<1e-5);
  }
  assert.ok(p.geometry.getAttribute('normal').array.every(Number.isFinite));p.geometry.dispose();
 }
 assert.equal(JSON.stringify(input),original);
});
test('oral and cranial scopes preserve source anatomy without embedding the independent pocket',()=>{
 const mouth=atlasContext('mouth'),face=atlasContext('face');
 assert.ok(face.length>mouth.length);
 assert.ok(mouth.some(p=>p.id==='palate'));assert.ok(mouth.some(p=>p.id==='tongue'));
 assert.equal(mouth.filter(p=>p.id==='gingiva').length,2);
 assert.ok(!mouth.some(p=>['plaque','lumen','pulp','epithelium'].includes(p.id)));
 assert.ok(ATLAS_FDI36.x>0);assert.ok(ATLAS_FDI36.toArray().every(Number.isFinite));
 for(const p of [...mouth,...face])p.geometry.dispose();
});
