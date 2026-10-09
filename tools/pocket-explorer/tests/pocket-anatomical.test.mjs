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
test('a three-tooth segment gives the gingiva somewhere natural to continue',()=>{
 const seg=asset.authored_features.segment;
 assert.deepEqual(seg.teeth,['FDI 35','FDI 36','FDI 37']);assert.equal(seg.selected,'FDI 36');
 const [x0,x1]=seg.ends;
 for(const [health,parts] of Object.entries(asset.cases)){
  assert.equal(parts.filter(p=>p.id==='enamel').length,3,health+' enamel of 35, 36 and 37');
  assert.ok(parts.some(p=>p.name.includes('FDI 35'))&&parts.some(p=>p.name.includes('FDI 37')),health);
  const gum=parts.filter(p=>p.id==='gingiva');assert.equal(gum.length,1,health+' one continuous gingiva');
  const g=gum[0].positions,bins=new Map();
  for(let i=0;i<g.length;i+=3){const x=g[i],z=g[i+2];if(x<x0-.01||x>x1+.01)continue;const key=Math.floor((x-x0)/(x1-x0)*80)+(z>0?':b':':l');bins.set(key,Math.min(bins.get(key)??Infinity,g[i+1]));}
  // Gingiva covers both the cheek and the tongue side along the whole segment.
  for(const side of [':b',':l'])for(let k=0;k<80;k++)assert.ok(bins.has(k+side),health+' gap in gingiva at bin '+k+side);
  // Its basal edge is a continuous rolled contour along each side, without the
  // Boolean teeth of the earlier voxel-repaired collar.
  for(const side of [':b',':l']){const low=Array.from({length:80},(_,k)=>bins.get(k+side));const worst=Math.max(...low.slice(1,-1).map((y,k)=>Math.abs(y-(low[k]+low[k+2])/2)));assert.ok(worst<0.25,health+side+' basal edge roughness '+worst.toFixed(3));}
 }
});
test('interdental papillae rise above the attached gingiva beside 36',()=>{
 const g=asset.cases.healthy.find(p=>p.id==='gingiva').positions;
 // Highest gingiva in a thin slice across the arch, between two |z| limits.
 const top=(x,z0,z1)=>{let y=-Infinity;for(let i=0;i<g.length;i+=3){const z=Math.abs(g[i+2]);if(Math.abs(g[i]-x)<.35&&z>=z0&&z<z1)y=Math.max(y,g[i+1]);}return y;};
 for(const x of [-5.4,5.2]){
  const papilla=top(x,0,1),attached=top(x,4.6,6.2);
  assert.ok(papilla>-1.5,'papilla at x='+x+' only reaches '+papilla.toFixed(2));
  assert.ok(papilla>attached+1,'papilla at x='+x+' ('+papilla.toFixed(2)+') not above attached gingiva ('+attached.toFixed(2)+')');
 }
});
test('pocket and junctional boundary retain separate shallow and deep extents',()=>{
 const ys=(health,id)=>asset.cases[health].find(p=>p.id===id).positions.filter((_,i)=>i%3===1);
 assert.ok(Math.min(...ys('healthy','lumen'))>Math.min(...ys('periodontitis','lumen'))+4);
 for(const health of ['healthy','periodontitis'])assert.ok(Math.min(...ys(health,'epithelium'))<Math.min(...ys(health,'lumen')));
 assert.equal(asset.authored_features.cusps.length,5);assert.equal(asset.authored_features.atlas_mesh_reuse,false);
});
