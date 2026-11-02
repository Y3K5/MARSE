import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {retainedGroups,explorationChoices,explorationOptions,sourceNetworkData,fieldLayout,foldRadius,displayLimits,supportedAcrossPairs,overviewPointGroups} from './exploration-model.js';
import {validateView} from './model.js';
import {sceneData,resultRows} from './subtractive-model.js';

const root=process.env.MARSE_RELATIONS_FIXTURE;
if(!root)throw Error('Set MARSE_RELATIONS_FIXTURE to invented software input.');
const data=JSON.parse(readFileSync(root+'/data/map.json'));

test('organism choices are metadata-driven; reference holds do not become folds or new review members',()=>{
 const before=JSON.stringify(data),groups=retainedGroups(data),choices=explorationChoices(data);
 assert(groups.length>0);assert(choices.every(c=>c.label&&!c.label.includes('undefined')));assert.equal(groups.reduce((n,g)=>n+g.review,0),data.subtractive.shortlist_keys.length);
 assert.equal(choices.find(c=>c.id==='review').records,3);assert.equal(choices.find(c=>c.id==='review').folds,2);
 assert.equal(choices.find(c=>c.id==='all').folds,2);
 const id=groups[0].id;assert.deepEqual(explorationOptions(data,'forms:'+id),{cohort:'forms',organism:id});assert.deepEqual(explorationOptions(data,'review'),{cohort:'shortlist',organism:'all'});assert.throws(()=>explorationOptions(data,'forms:unknown-homolog'));
 assert.equal(JSON.stringify(data),before);
});

test('the complete source atlas survives filtering out the retained protein field',()=>{
 const source=sourceNetworkData(data),ids=new Set([...source.proteins,...source.contexts].map(p=>p.id));
 assert(source.relationships.every(r=>ids.has(r.a)&&ids.has(r.b)));
 assert.deepEqual(source.relationships,data.relationships);assert(source.proteins.every(p=>data.subtractive.source_context_ids.includes(p.id)));
 assert(!source.proteins.some(p=>p.subtractive_record));
});

function inventedPopulation(){
 const species=['alpha','beta','gamma'],proteins=[],shape_traces={};
 for(const [s,count]of species.map((s,i)=>[s,[37,61,83][i]]))for(let i=0;i<count;i++){const id=s+'-'+i;proteins.push({id,species:s,name:id,subtractive_record:{primary_group:'group-'+i%4}});shape_traces[id]={};}
 return {species:species.map(id=>({id})),proteins,contexts:[],shape_traces};
}
test('dense fields keep positive envelope clearance within and across groups at every control extreme',()=>{
 const sample=inventedPopulation(),before=JSON.stringify(sample),centers={alpha:[-30,-30,0],beta:[30,0,20],gamma:[-30,30,0]};
 for(const mode of ['cluster','layers','orbit'])for(const size of [1,displayLimits.foldDefault,displayLimits.foldMax])for(const gap of [displayLimits.spacingMin,displayLimits.spacingMax]){
  const layout=fieldLayout(sample,gap,mode,size,centers),nodes=Object.values(layout);assert.equal(nodes.length,sample.proteins.length);
  for(const node of nodes)assert(node.position.every(Number.isFinite));
  for(let i=0;i<nodes.length;i++)for(let j=0;j<i;j++){
   const distance=Math.hypot(...nodes[i].position.map((v,k)=>v-nodes[j].position[k]));
   assert(distance>2*foldRadius(size)+4,mode+' display envelopes overlap');
  }
  assert(new Set(nodes.filter(n=>n.species==='alpha').map(n=>n.position[1])).size>1,'Field lost depth variation');
 }
 assert.equal(JSON.stringify(sample),before);
});

test('a focused reference is centered without empty organism anchors widening the field',()=>{
 const sample=inventedPopulation();sample.proteins=sample.proteins.filter(p=>p.species==='alpha');
 const layout=fieldLayout(sample,1.2,'layers',3.2,{alpha:[-45,0,0],beta:[45,0,0],gamma:[0,45,0]});
 assert(Object.values(layout).every(n=>n.groupCenter.every(v=>v===0)));
 assert.equal(Object.keys(layout).length,37);
});

test('saved views accept the larger display limits and safe optional overview camera; old views stay valid',()=>{
 const view={schema:'marse-lattice-view',version:1,sourceSHA256:'invented-input',preset:'microbes',context:{kind:'network'},spacing:3,explode:0,quality:'high',grid:true,labels:true,cutaway:false,camera:{position:[20,10,20],target:[0,0,0]},networkDisplay:{scope:'all',candidates:false,bodies:false,unlinked:true,layout:'layers',foldSize:4.5,overviewCamera:{position:[30,10,40],target:[0,0,0]}}};
 assert.equal(validateView(view,{sha256:'invented-input',map:data}),view);
 for(const change of [v=>v.spacing=3.1,v=>v.networkDisplay.foldSize=4.6,v=>v.networkDisplay.overviewCamera.position[0]=Infinity,v=>v.networkDisplay.overviewCamera.target=[30,10,40]]){
  const bad=structuredClone(view);change(bad);assert.throws(()=>validateView(bad,{sha256:'invented-input',map:data}));
 }
 delete view.networkDisplay.overviewCamera;view.spacing=1;view.networkDisplay.foldSize=1.8;assert.equal(validateView(view,{sha256:'invented-input',map:data}),view);
});


test('combined choice retains original forms and literal source endpoints without promoting reference links into interactions',()=>{
 const before=JSON.stringify(data),options=explorationOptions(data,'together'),choice=explorationChoices(data).find(c=>c.id==='together');
 assert.deepEqual(options,{cohort:'forms',organism:'all',sourceContext:true});
 const visible=sceneData(data,options),ids=new Set([...visible.proteins,...visible.contexts].map(p=>p.id));
 for(const p of resultRows(data,{cohort:'forms'}))assert(ids.has(p.id));
 for(const id of data.subtractive.source_context_ids)assert(ids.has(id));
 for(const p of data.proteins.filter(p=>data.subtractive.shortlist_keys.includes(p.subtractive_record?.key)))assert(ids.has(p.id),'A review coordinate gap disappeared from the combined scene');
 assert.deepEqual(visible.relationships,data.relationships);
 assert.equal(choice.folds,2);assert.equal(choice.sourceFolds,3);
 assert(visible.relationships.every(r=>data.relationships.includes(r)));
 assert.equal(JSON.stringify(data),before);
});

test('only literal named-protein direct source pairs across organisms get supported emphasis',()=>{
 const sample={proteins:[{id:'a',species:'alpha'},{id:'b',species:'beta'},{id:'same',species:'alpha'},{id:'unmapped-homolog',species:'beta'}],contexts:[{id:'cell',species:'beta'},{id:'class',species:'beta'}],relationships:[]};
 const types=['direct assay','direct biochemical','biochemical interaction','indirect / phenotype','direct cell-binding / receptor unresolved','community phenotype','binding preparation / context dependent','candidate / unresolved identity','structural and biochemical'];
 sample.relationships=types.map((level,i)=>({id:'r'+i,a:'a',b:'b',level}));
 sample.relationships.push({id:'cell-assay',a:'a',b:'cell',level:'direct assay'},{id:'class-cut',a:'a',b:'class',level:'direct biochemical'},{id:'within',a:'a',b:'same',level:'direct assay'});
 sample.subtractive={reference_links:[{a:'unmapped-homolog',b:'b',kind:'exact_reference_identity'}]};
 const before=JSON.stringify(sample),supported=supportedAcrossPairs(sample);
 assert.deepEqual(supported.map(r=>r.id),['r0','r1']);
 assert(supported.every(r=>r.a==='a'&&r.b==='b')); // No counterpart or homolog substitution.
 assert.equal(JSON.stringify(sample),before);
});


test('packed combined fields keep every display envelope separate in all layouts and control extremes',()=>{
 const sample=inventedPopulation();sample.display_field='combined';
 for(let i=0;i<sample.proteins.length;i+=9)delete sample.proteins[i].subtractive_record;
 const before=JSON.stringify(sample),centers={alpha:[-30,-30,0],beta:[30,0,20],gamma:[-30,30,0]};
 for(const mode of ['cluster','layers','orbit'])for(const size of [displayLimits.foldMin,displayLimits.foldMax])for(const gap of [displayLimits.spacingMin,displayLimits.spacingMax]){
  const nodes=Object.values(fieldLayout(sample,gap,mode,size,centers));assert.equal(nodes.length,181);
  for(let i=0;i<nodes.length;i++)for(let j=0;j<i;j++)assert(Math.hypot(...nodes[i].position.map((v,k)=>v-nodes[j].position[k]))>2*foldRadius(size)+4);
 }
 assert.equal(JSON.stringify(sample),before);
});

test('the source overview draws an unassigned context endpoint without presenting it as another microbe',()=>{
 const sample={proteins:[{id:'known',name:'Known source',species:'alpha'}],contexts:[{id:'unresolved',name:'Unresolved endpoint',species:'unknown',kind:'context'}],shape_traces:{},relationships:[{id:'context-link',a:'known',b:'unresolved'}],subtractive:{source_context_ids:['known']}};
 const before=JSON.stringify(sample),source=sourceNetworkData(sample),layout=fieldLayout(source,.7,'cluster',1,{alpha:[0,0,0],unknown:[200,0,0]}),groups=overviewPointGroups(layout);
 assert.deepEqual(groups.flatMap(g=>g.nodes.map(n=>n.id)).sort(),['known','unresolved']);assert.equal(groups.find(g=>g.species==='unknown').nodes[0].kind,'context');
 assert.equal(source.relationships.length,1);assert.equal(JSON.stringify(sample),before);
});
