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
test('basal gingival edge is a continuous rolled contour without Boolean teeth',()=>{
 // Lowest gingival vertex per azimuth traces the basal edge on the cortical plate.
 // The earlier voxel-repaired collar measured 1.15 / 0.80 model units here.
 for(const [health,parts] of Object.entries(asset.cases)){
  const g=parts.find(p=>p.id==='gingiva'),bins=180,low=Array(bins).fill(Infinity);
  for(let i=0;i<g.positions.length;i+=3){const b=Math.floor((Math.atan2(g.positions[i+2],g.positions[i])+Math.PI)/(2*Math.PI)*bins)%bins;low[b]=Math.min(low[b],g.positions[i+1]);}
  assert.ok(low.every(Number.isFinite),health+' gingiva leaves an azimuth uncovered');
  const worst=Math.max(...low.map((y,i)=>Math.abs(y-(low[(i+bins-1)%bins]+low[(i+1)%bins])/2)));
  assert.ok(worst<0.25,health+' basal edge roughness '+worst.toFixed(3));
 }
});
test('every exported triangle edge is shared by exactly two triangles',()=>{
 for(const [health,parts] of Object.entries(asset.cases))for(const part of parts){
  const uses=new Map();for(let i=0;i<part.indices.length;i+=3)for(let k=0;k<3;k++){const a=part.indices[i+k],b=part.indices[i+(k+1)%3],key=a<b?a*65536+b:b*65536+a;uses.set(key,(uses.get(key)||0)+1);}
  assert.equal([...uses.values()].filter(c=>c!==2).length,0,health+' '+part.name);
 }
});
test('plaque compartments follow the tooth above and below the margin without species positions',()=>{
 for(const [health,parts] of Object.entries(asset.cases)){
  const supra=parts.filter(p=>p.id==='supragingival'),sub=parts.filter(p=>p.id==='plaque');
  assert.equal(supra.length,1,health);assert.ok(sub.length>=2,health+' needs the site ribbon and the circumferential film');
  // The circumferential film covers both the cheek (+z) and tongue (-z) sides.
  const film=sub.find(p=>p.name.includes('Circumferential')),z=film.positions.filter((_,i)=>i%3===2);
  assert.ok(Math.max(...z)>3&&Math.min(...z)<-3,health);
  for(const p of [...supra,...sub])assert.ok(volume(p)>0.05&&volume(p)<5,health+' '+p.name);
 }
});
test('pocket and junctional boundary retain separate shallow and deep extents',()=>{
 const ys=(health,id)=>asset.cases[health].find(p=>p.id===id).positions.filter((_,i)=>i%3===1);
 assert.ok(Math.min(...ys('healthy','lumen'))>Math.min(...ys('periodontitis','lumen'))+4);
 for(const health of ['healthy','periodontitis'])assert.ok(Math.min(...ys(health,'epithelium'))<Math.min(...ys(health,'lumen')));
 assert.equal(asset.authored_features.cusps.length,5);assert.equal(asset.authored_features.atlas_mesh_reuse,false);
});
