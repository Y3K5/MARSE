import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pairLimits,validatePairFocus,pairSceneData,pairLayout} from './pair-model.js';
import {foldRadius,displayLimits} from './exploration-model.js';
import {validateView} from './model.js';
const root=process.env.MARSE_RELATIONS_FIXTURE;if(!root)throw Error('Set MARSE_RELATIONS_FIXTURE to invented software inputs.');
const data=JSON.parse(readFileSync(root+'/data/map.json')),state={record:'pair-1',magnification:pairLimits.default};

test('large pair keeps the literal supported source endpoints; missing coordinates stay missing',()=>{
 const before=JSON.stringify(data),record=validatePairFocus(data,state),scene=pairSceneData(data,state);
 assert.deepEqual(scene.proteins.map(p=>p.id).sort(),[record.a,record.b].sort());assert.deepEqual(scene.relationships,[record]);assert.equal(scene.contexts.length,0);
 assert(!scene.proteins.some(p=>p.subtractive_record),'A reference connector was promoted into a source endpoint');
 assert.equal(JSON.stringify(data),before);
 const missing=structuredClone(data);delete missing.shape_traces[record.a];const missingBefore=JSON.stringify(missing),layout=pairLayout(missing,state,3.2,1.5);
 assert.equal(layout[record.a].hasShape,false);assert.equal(layout[record.b].hasShape,true);assert.equal(JSON.stringify(missing),missingBefore);
 for(const bad of [{...state,record:'pair-2'},{...state,record:'unknown'},{...state,magnification:0},{...state,magnification:4.1},{...state,magnification:'3'},{...state,magnification:Infinity},{...state,endpoint:'guessed'},{...state,path:'/not-allowed'}])assert.throws(()=>validatePairFocus(data,bad));
});

test('both displayed folds enlarge with positive envelope clearance at every control extreme',()=>{
 const before=JSON.stringify(data);
 for(const foldSize of [displayLimits.foldMin,displayLimits.foldDefault,displayLimits.foldMax])for(const spacing of [displayLimits.spacingMin,displayLimits.spacingMax])for(const magnification of [pairLimits.min,pairLimits.default,pairLimits.max]){
  const layout=pairLayout(data,{...state,magnification},foldSize,spacing),nodes=Object.values(layout);
  assert.equal(nodes.length,2);assert(nodes.every(n=>n.displayRadius===foldRadius(foldSize)*magnification));
  const distance=Math.hypot(...nodes[0].position.map((v,k)=>v-nodes[1].position[k]));assert(distance>nodes[0].displayRadius+nodes[1].displayRadius+4);
  assert.equal(nodes[0].position[0],-nodes[1].position[0]);assert(nodes.every(n=>n.position[1]===0&&n.position[2]===0));
 }
 for(const args of [[0,1],[5,1],[2,.1],[2,4]])assert.throws(()=>pairLayout(data,state,...args));
 assert.equal(JSON.stringify(data),before);
});

test('saved pair choice is source-pinned and bounded while legacy views remain valid',()=>{
 const view={schema:'marse-lattice-view',version:1,sourceSHA256:'invented-source',preset:'microbes',context:{kind:'network'},spacing:1.5,explode:0,quality:'high',grid:true,labels:true,cutaway:false,camera:{position:[0,10,100],target:[0,0,0]},networkDisplay:{scope:'across',candidates:false,bodies:false,unlinked:true,layout:'layers',foldSize:3.2,pairFocus:state,subtractive:{cohort:'forms',organism:'all',group:'all',query:'',page:0,sourceContext:true,referenceLinks:true}}};
 const pinned={sha256:'invented-source',map:data};assert.equal(validateView(view,pinned),view);
 for(const edit of [v=>v.networkDisplay.pairFocus.record='pair-2',v=>v.networkDisplay.pairFocus.magnification=5,v=>v.networkDisplay.subtractive.sourceContext=false,v=>v.networkDisplay.subtractive.cohort='shortlist',v=>v.context.kind='overview',v=>v.sourceSHA256='stale']){const bad=structuredClone(view);edit(bad);assert.throws(()=>validateView(bad,pinned));}
 const legacy=structuredClone(view);delete legacy.networkDisplay.pairFocus;assert.equal(validateView(legacy,pinned),legacy);
});
