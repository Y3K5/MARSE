import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {navigationItems,navigationCursor,proteinNavigationItem,validateNavigationProtein} from './navigation-model.js';
import {validateView} from './model.js';
const root=process.env.MARSE_RELATIONS_FIXTURE;if(!root)throw Error('Set MARSE_RELATIONS_FIXTURE to invented software inputs.');
const data=JSON.parse(readFileSync(root+'/data/map.json'));

test('navigation follows original filtered metadata/order and identifies coordinate holds without adopting targets',()=>{
 const before=JSON.stringify(data),items=navigationItems(data,{cohort:'shortlist',organism:'all',sourceContext:false});
 assert.equal(items.length,3);assert(items.every(i=>i.name&&i.identity&&i.scope==='Retained result'));assert.equal(items.filter(i=>!i.hasShape).length,1);assert.match(items.find(i=>!i.hasShape).status,/no fold assigned/);
 assert.equal(navigationCursor(items,null).next,items[0].id);assert.equal(navigationCursor(items,items[0].id).previous,null);assert.equal(navigationCursor(items,items[0].id).next,items[1].id);assert.equal(navigationCursor(items,items.at(-1).id).next,null);assert.deepEqual(navigationCursor([],null),{index:-1,total:0,previous:null,next:null});
 assert.equal(JSON.stringify(data),before);
});

test('pair navigation retains literal source names and source scope instead of projecting onto retained references',()=>{
 const before=JSON.stringify(data),items=navigationItems(data,{},true);assert.equal(items.length,1);assert.equal(items[0].id,'pair-1');assert.equal(items[0].scope,'Published-source relationship');assert(items[0].name.includes(data.proteins.find(p=>p.id==='source-a').name));assert(items[0].name.includes(data.proteins.find(p=>p.id==='source-b').name));
 assert.equal(proteinNavigationItem(data,'source-a').scope,'Source endpoint');assert.equal(proteinNavigationItem(data,'not-present'),null);assert.equal(JSON.stringify(data),before);
});

test('selection must be present in the actual reference/filter/page; a reference link is not visibility or identity substitution',()=>{
 const options={cohort:'shortlist',organism:'all',sourceContext:false},items=navigationItems(data,options);assert.equal(validateNavigationProtein(data,options,items[0].id),items[0].id);
 assert.throws(()=>validateNavigationProtein(data,options,'source-a'));assert.throws(()=>validateNavigationProtein(data,{...options,query:'no-such-protein'},items[0].id));assert.throws(()=>validateNavigationProtein(data,{...options,organism:'unrelated-reference'},items[0].id));assert.throws(()=>validateNavigationProtein(data,options,'unknown'));
 assert.equal(validateNavigationProtein(data,{...options,sourceContext:true},'source-a'),'source-a');
});

test('source-pinned saved current names refuse conflicting frames/filter scopes and preserve legacy saves',()=>{
 const options={cohort:'shortlist',organism:'all',group:'all',query:'',page:0,sourceContext:false,referenceLinks:true},id=navigationItems(data,options)[0].id,view={schema:'marse-lattice-view',version:1,sourceSHA256:'invented-source',preset:'microbes',context:{kind:'network'},spacing:1.5,explode:0,quality:'balanced',grid:true,labels:true,cutaway:false,camera:{position:[0,10,100],target:[0,0,0]},networkDisplay:{scope:'all',candidates:false,bodies:false,unlinked:true,layout:'layers',foldSize:3.2,subtractive:options,navigationProtein:id}},pinned={sha256:'invented-source',map:data};
 assert.equal(validateView(view,pinned),view);
 for(const edit of [v=>v.networkDisplay.navigationProtein='source-a',v=>v.networkDisplay.navigationProtein='unknown',v=>v.networkDisplay.subtractive.query='no-such-record',v=>v.context.kind='overview',v=>v.sourceSHA256='stale',v=>{v.networkDisplay.subtractive.cohort='forms';v.networkDisplay.subtractive.sourceContext=true;v.networkDisplay.pairFocus={record:'pair-1',magnification:2.4};}]){const bad=structuredClone(view);edit(bad);assert.throws(()=>validateView(bad,pinned));}
 const legacy=structuredClone(view);delete legacy.networkDisplay.navigationProtein;assert.equal(validateView(legacy,pinned),legacy);
});
