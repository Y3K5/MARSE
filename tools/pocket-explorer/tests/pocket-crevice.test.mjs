import {test} from 'node:test';
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {readFile} from 'node:fs/promises';
import {buildCrevice,microbialObject,microbialPlacements,creviceGap,creviceLandmarks,creviceCamera,CREVICE_CLAIM} from '../dist/pocket-core/pocket-crevice.mjs';
import {TAXA,PRESETS} from '../dist/pocket-core/pocket-science.mjs';
const hash=g=>createHash('sha256').update(Buffer.from(g.getAttribute('position').array.buffer)).digest('hex');
test('all nine objects are deterministic, distinct, finite and source identified',()=>{
 const hashes=[];
 for(const t of TAXA){const a=microbialObject(t.id),b=microbialObject(t.id);assert.equal(hash(a),hash(b));hashes.push(hash(a));for(const key of ['position','normal'])assert.ok(a.getAttribute(key).array.every(Number.isFinite));a.dispose();b.dispose();}
 assert.equal(new Set(hashes).size,9);assert.throws(()=>microbialObject('unresolved'),/Unknown microbial/);
 for(const health of ['healthy','periodontitis']){const parts=buildCrevice(health),cells=parts.filter(p=>p.taxon);assert.deepEqual(cells.map(c=>c.taxon),TAXA.map(t=>t.id));for(const c of cells)assert.match(c.source,/^https:\/\//);for(const p of parts){assert.ok(p.geometry.getAttribute('normal').array.every(Number.isFinite));p.geometry.dispose();}}
});
test('entrance stays narrower than the authored interior and landmarks remain separate',()=>{
 const original=JSON.stringify(PRESETS);
 for(const health of ['healthy','periodontitis']){const {margin,attachment,crest}=creviceLandmarks(health);assert.ok(margin>attachment&&attachment>crest);assert.ok(creviceGap(margin,health)<creviceGap((margin+attachment)/2,health));assert.ok(creviceGap(attachment,health)<creviceGap((margin+attachment)/2,health));for(const p of microbialPlacements(health))assert.ok(p.position[1]<margin&&p.position[1]>attachment);}
 assert.equal(JSON.stringify(PRESETS),original);
});
test('crevice cameras are distinct and the isolated form focus follows selected identity',()=>{
 const poses=['section','entrance','wall'].map(v=>creviceCamera(v,'periodontitis'));assert.equal(new Set(poses.map(p=>JSON.stringify(p))).size,3);
 for(const p of poses)assert.ok([...p.position,...p.target].every(Number.isFinite));
 assert.deepEqual(creviceCamera('wall','periodontitis','td').target,microbialPlacements('periodontitis').find(p=>p.id==='td').position);
});
test('3D examples are independent from immutable saved endpoints and deny measured placement',async()=>{
 const path=new URL('../src/pocket-input.json',import.meta.url),before=await readFile(path);for(const health of ['healthy','periodontitis'])for(const p of buildCrevice(health))p.geometry.dispose();assert.deepEqual(await readFile(path),before);assert.match(CREVICE_CLAIM,/not measured colonization/);assert.match(CREVICE_CLAIM,/sizes exaggerated independently/);
});
test('every authored tissue loft has closed, consistently oriented section boundaries',()=>{
 for(const health of ['healthy','periodontitis'])for(const part of buildCrevice(health)){
  const g=part.geometry;if(part.taxon){g.dispose();continue;}
  const a=g.getAttribute('position'),ix=g.index.array,edges=new Map();let volume=0;
  const key=i=>[a.getX(i),a.getY(i),a.getZ(i)].map(n=>Math.round(n*1e5)).join(',');
  for(let j=0;j<ix.length;j+=3){
   const p=Array.from({length:3},(_,k)=>[a.getX(ix[j+k]),a.getY(ix[j+k]),a.getZ(ix[j+k])]);
   const u=p[1].map((n,k)=>n-p[0][k]),v=p[2].map((n,k)=>n-p[0][k]),cross=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];
   assert.ok(Math.hypot(...cross)>1e-9,health+' '+part.id+' degenerate triangle');
   volume+=p[0][0]*(p[1][1]*p[2][2]-p[1][2]*p[2][1])+p[0][1]*(p[1][2]*p[2][0]-p[1][0]*p[2][2])+p[0][2]*(p[1][0]*p[2][1]-p[1][1]*p[2][0]);
   for(let k=0;k<3;k++){const x=key(ix[j+k]),y=key(ix[j+(k+1)%3]),e=x<y?x+'|'+y:y+'|'+x,count=edges.get(e)||[0,0];count[0]++;count[1]+=x<y?1:-1;edges.set(e,count);}
  }
  assert.ok([...edges.values()].every(v=>v[0]===2&&v[1]===0),health+' '+part.id+' open or reversed seam');assert.ok(volume>0);
  assert.deepEqual(g.groups.map(p=>p.materialIndex),[0,1]);g.dispose();
 }
});
