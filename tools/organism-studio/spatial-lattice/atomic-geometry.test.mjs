// Optional renderer checks use the already-installed local Three assets, never remote code.
import assert from 'node:assert/strict';
import test from 'node:test';
import {register} from 'node:module';
const three=new URL('../vendor/three.module.js',import.meta.url).href;
register('data:text/javascript,'+encodeURIComponent(`export async function resolve(s,c,next){if(s==='three')return {url:${JSON.stringify(three)},shortCircuit:true};return next(s,c);}`));
const T=await import('three'),{AtomicLens}=await import('./atomic-lens.js');

test('repeated atom selection retains exact local geometry and disposes every replaced instance buffer',()=>{
 const payload={atoms:[0,1.5,3].map((x,i)=>({i,xyz:[x,0,0],chain:'A',element:'C'})),bonds:[{a:0,b:1,order:'sing',kind:'component_template'},{a:1,b:2,order:'doub',kind:'component_template'}]},before=JSON.stringify(payload),lens=Object.create(AtomicLens.prototype);
 Object.assign(lens,{payload,selected:0,chain:'all',environment:true,scene:new T.Scene(),canvas:{dataset:{atomicState:'old'}},onAtom:()=>{}});lens.draw();
 let disposed=0;for(let n=0;n<40;n++){
  lens.root.traverse(o=>{if(o.isInstancedMesh)o.addEventListener('dispose',()=>disposed++);});
  lens.select(n%3);assert.equal(lens.scene.children.length,1);assert.equal(lens.visible.length,3);assert.equal(lens.drawnBonds,2);
  for(const [i,index]of lens.visible.entries()){
   const matrix=new T.Matrix4();lens.spheres.getMatrixAt(i,matrix);const pos=new T.Vector3().setFromMatrixPosition(matrix),expected=payload.atoms[index].xyz.map((v,k)=>v-payload.atoms[lens.selected].xyz[k]);assert.deepEqual(pos.toArray(),expected);
  }
 }
 assert.equal(disposed,80);assert.equal(JSON.stringify(payload),before);lens.clear();assert.equal(lens.payload,null);assert.equal(lens.scene.children.length,0);assert.equal(lens.canvas.dataset.atomicState,undefined);
});
