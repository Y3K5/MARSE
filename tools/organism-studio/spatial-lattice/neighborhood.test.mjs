import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {proteinNeighborhood,neighborhoodScene,neighborhoodLayout,neighborhoodSites,validateNeighborhood,validateNeighborhoodSelection,neighborhoodReturnView} from './neighborhood-model.js';
import {loadRelationPacket} from './interaction-model.js';
import {sourceNetworkData,displayLimits} from './exploration-model.js';
import {validateView} from './model.js';
const root=process.env.MARSE_RELATIONS_FIXTURE;if(!root)throw Error('Set MARSE_RELATIONS_FIXTURE to invented software inputs.');
const read=n=>readFileSync(root+'/'+n),json=n=>JSON.parse(read(n)),hash=b=>createHash('sha256').update(b).digest('hex'),data=json('data/map.json'),ledger=json('residues.json');
await loadRelationPacket(data,ledger,hash(read('data/map.json')),hash(read('residues.json')),async n=>new Response(read(n)));
test('retained reference, exact source endpoint and literal supported partner stay separate; whole atlas is unchanged',()=>{
 const before=JSON.stringify(data),state={protein:'review-a'},n=proteinNeighborhood(data,state),scene=neighborhoodScene(data,state),all=sourceNetworkData(data);
 assert.deepEqual(n.sourceIds,['source-a']);assert.deepEqual(n.records.map(r=>r.id),['pair-1']);assert.deepEqual(n.partners,['source-b']);assert.deepEqual(scene.proteins.map(p=>p.id),['review-a','source-a','source-b']);assert(!scene.relationships.some(r=>r.a==='review-a'||r.b==='review-a'));assert.equal(n.references.length,1);assert.equal(scene.contexts.length,0);assert.equal(all.relationships.length,data.relationships.length);assert.equal(JSON.stringify(data),before);
});
test('unmatched sequence, species and label-only holds cannot gain a partner; coordinate gaps stay markers',()=>{
 for(const edit of [d=>d.proteins.find(p=>p.id==='review-a').subtractive_record.sequence_sha256='wrong',d=>d.proteins.find(p=>p.id==='source-a').species='beta',d=>d.subtractive.reference_links[0].kind='homolog',d=>d.subtractive.reference_links[0].sequence_sha256='']){const d=structuredClone(data);edit(d);const n=proteinNeighborhood(d,{protein:'review-a'});assert.equal(n.records.length,0);assert.deepEqual(n.proteins.map(p=>p.id),['review-a']);}
 const held=proteinNeighborhood(data,{protein:'review-hold'});assert.equal(held.records.length,0);assert.equal(held.holds.length,1);assert.equal(held.references.length,0);
 const gap=neighborhoodLayout(data,{protein:'review-gap'},3.2,1.5);assert.equal(Object.keys(gap).length,1);assert.equal(gap['review-gap'].hasShape,false);
 for(const state of [null,false,[],{protein:'missing'},{protein:data.contexts[0].id},{protein:'source-a',path:'/not-accepted'}])assert.throws(()=>validateNeighborhood(data,state));
});
test('cell/class, candidate, phenotype and association records stay outside supported protein edges',()=>{
 const d=structuredClone(data);for(const [i,level] of ['biochemical interaction','community phenotype','indirect / phenotype','candidate / unresolved identity'].entries())d.relationships.push({id:'other-'+i,a:'source-a',b:'source-c',level});
 d.relationships.push({id:'unknown-cell',a:'source-a',b:d.contexts[0].id,level:'direct assay'});const n=proteinNeighborhood(d,{protein:'source-a'});assert.deepEqual(n.records.map(r=>r.id),['pair-1']);assert(n.otherRecords.some(r=>r.id==='unknown-cell'));assert(!n.otherRecords.some(r=>r.level==='candidate / unresolved identity'));assert(!n.proteins.some(p=>p.id==='source-c'));
});
test('site membership comes from the qualified reader on literal source folds, never on retained seed or unsupported partner',()=>{
 const before=JSON.stringify(data),sites=neighborhoodSites(data,ledger,{protein:'review-a'});assert.equal(sites.length,1);assert.equal(sites[0].endpoint,'source-a');assert.equal(sites[0].region.category,'assayed_construct');assert.equal(sites[0].region.correspondence.length,2);assert(!sites.some(s=>s.endpoint==='review-a'||s.endpoint==='source-b'));assert.deepEqual(neighborhoodSites(data,ledger,{protein:'review-hold'}),[]);assert.equal(JSON.stringify(data),before);
});
test('neighborhood folds remain separated in true XYZ across every display extreme and all layouts',()=>{
 const d=structuredClone(data);for(let i=0;i<35;i++){d.proteins.push({id:'extra-'+i,name:'Invented neighbor',species:'beta'});d.relationships.push({id:'invented-'+i,a:'source-a',b:'extra-'+i,level:'direct assay'});}
 const before=JSON.stringify(d);for(const mode of ['cluster','layers','orbit'])for(const size of [1,displayLimits.foldDefault,displayLimits.foldMax])for(const spacing of [displayLimits.spacingMin,displayLimits.spacingMax]){const nodes=Object.values(neighborhoodLayout(d,{protein:'review-a'},size,spacing,mode));assert.equal(nodes.length,38);assert.deepEqual(nodes[0].position,[0,0,0]);assert(new Set(nodes.map(n=>n.position[1])).size>1);assert(new Set(nodes.map(n=>n.position[2])).size>1);for(let i=0;i<nodes.length;i++)for(let j=0;j<i;j++)assert(Math.hypot(...nodes[i].position.map((v,k)=>v-nodes[j].position[k]))>nodes[i].displayRadius+nodes[j].displayRadius+4);}
 assert.equal(JSON.stringify(d),before);
});
test('saved neighborhood is source-pinned, rejects unrelated selection/pair/frame state and preserves old saves',()=>{
 const v={schema:'marse-lattice-view',version:1,sourceSHA256:'invented',preset:'microbes',context:{kind:'network'},spacing:1.5,explode:0,quality:'high',grid:true,labels:true,cutaway:false,camera:{position:[10,10,100],target:[0,0,0]},networkDisplay:{scope:'all',candidates:false,bodies:false,unlinked:true,neighborhood:{protein:'review-a'},navigationProtein:'source-b'}},pin={sha256:'invented',map:data};assert.equal(validateView(v,pin),v);assert.equal(validateNeighborhoodSelection(data,v.networkDisplay.neighborhood,'source-a'),'source-a');
 for(const edit of [s=>s.sourceSHA256='stale',s=>s.networkDisplay.neighborhood.protein='unknown',s=>s.networkDisplay.neighborhood=null,s=>s.networkDisplay.navigationProtein='source-c',s=>s.context.kind='overview',s=>s.context={kind:'protein',id:'review-a'},s=>s.networkDisplay.pairFocus={record:'pair-1',magnification:2.4}]){const bad=structuredClone(v);edit(bad);assert.throws(()=>validateView(bad,pin));}
 const old=structuredClone(v);delete old.networkDisplay.neighborhood;delete old.networkDisplay.navigationProtein;assert.equal(validateView(old,pin),old);
});

test('return-to-field camera is captured only from the molecular display field, never from source-coordinate or designed frames',()=>{
 const camera={position:[10,20,100],target:[1,2,3]},before=JSON.stringify(camera),saved=neighborhoodReturnView('microbes',{kind:'network'},camera);assert.deepEqual(saved,camera);assert.notEqual(saved.position,camera.position);
 for(const kind of ['protein','organism','overview'])assert.equal(neighborhoodReturnView('microbes',{kind},camera),null);assert.equal(neighborhoodReturnView('desk',{kind:'network'},camera),null);assert.equal(JSON.stringify(camera),before);
});
