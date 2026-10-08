import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
const asset=JSON.parse(readFileSync(new URL('../assets/pocket-organic/anatomy.json',import.meta.url)));
const volume=p=>{let sum=0;const v=p.positions;for(let i=0;i<p.indices.length;i+=3){const [a,b,c]=p.indices.slice(i,i+3).map(n=>v.slice(n*3,n*3+3));sum+=(a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0]))/6;}return sum;};
test('supporting compartments are substantial closed volumes, with two distinct bone regions',()=>{
 for(const [health,parts] of Object.entries(asset.cases)){
  for(const part of parts){assert.ok(part.indices.length>=36,health+' '+part.name+' empty');assert.ok(volume(part)>0.001,health+' '+part.name+' degenerate');}
  const bones=parts.filter(p=>p.id==='bone');assert.equal(bones.length,2);
  assert.ok(bones.some(p=>p.name.includes('cortical')));assert.ok(bones.some(p=>p.name.includes('Trabecular')));
  for(const bone of bones)assert.ok(volume(bone)>100,health+' '+bone.name+' collapsed');
  for(const id of ['pdl','connective','gingiva'])assert.ok(volume(parts.find(p=>p.id===id))>5,health+' '+id+' collapsed');
 }
});
test('pocket and junctional boundary retain separate shallow and deep extents',()=>{
 const ys=(health,id)=>asset.cases[health].find(p=>p.id===id).positions.filter((_,i)=>i%3===1);
 assert.ok(Math.min(...ys('healthy','lumen'))>Math.min(...ys('periodontitis','lumen'))+4);
 for(const health of ['healthy','periodontitis'])assert.ok(Math.min(...ys(health,'epithelium'))<Math.min(...ys(health,'lumen')));
 assert.equal(asset.authored_features.cusps.length,5);assert.equal(asset.authored_features.atlas_mesh_reuse,false);
});
