import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import test from 'node:test';
import {joinSiteResidues,siteResidueEntries,siteResiduePage,siteLabelCandidates,residueEvidenceClaim} from './site-residue-model.js';
import {loadRelationPacket,endpointLocation} from './interaction-model.js';
const root=process.env.MARSE_RELATIONS_FIXTURE;if(!root)throw Error('Set MARSE_RELATIONS_FIXTURE to invented software inputs.');
const read=n=>readFileSync(root+'/'+n),json=n=>JSON.parse(read(n)),hash=b=>createHash('sha256').update(b).digest('hex'),data=json('data/map.json'),ledger=json('residues.json');
await loadRelationPacket(data,ledger,hash(read('data/map.json')),hash(read('residues.json')),async n=>new Response(read(n)));

test('labels keep independent paper, canonical and author/chain/insertion identities; numbers alone never join',()=>{
 const rows=[{index:0,chain:'A',number:368,insertion:'',name:'LYS',point:[368,1,2,3]},{index:1,chain:'B',number:368,insertion:'',name:'LYS',point:[368,3,2,1]},{index:2,chain:'A',number:368,insertion:'A',name:'THR',point:[368,4,5,6]}],correspondence=[{paper:368,canonical:392,author:368,chain:'A',insertion:'',name:'LYS'},{paper:369,canonical:393,author:368,chain:'A',insertion:'A',name:'THR'},{paper:100,canonical:200,author:368,chain:'B',insertion:'',name:'LYS'}],before=JSON.stringify({rows,correspondence}),joined=joinSiteResidues(rows,correspondence);
 assert.deepEqual(joined.map(r=>r.index),[0,2,1]);assert.equal(joined[0].canonical,392);assert.equal(joined[0].authorLabel,'A:368');assert.equal(joined[1].authorLabel,'A:368A');assert.equal(joined[2].authorLabel,'B:368');
 for(const edit of [c=>c[0].name='GLY',c=>c[0].chain='C',c=>c[0].insertion='B',c=>c[0].author=369,c=>c[0].canonical=null,c=>c.push({...c[0]})]){const bad=structuredClone(correspondence);edit(bad);assert.throws(()=>joinSiteResidues(rows,bad));}
 assert.throws(()=>joinSiteResidues(rows.slice(1),correspondence));assert.throws(()=>joinSiteResidues([...rows,rows[0]],correspondence));const missing=structuredClone(rows);missing[0].point=null;assert.throws(()=>joinSiteResidues(missing,correspondence));
 assert.equal(JSON.stringify({rows,correspondence}),before);
});

test('site entries consume the existing qualified reader and refuse stale selection; unsupported partner stays empty',()=>{
 const before=JSON.stringify(data),loc=endpointLocation(data,ledger,'pair-1','source-a'),region=loc.regions[0],entries=siteResidueEntries(data,ledger,loc,region);
 assert.equal(entries.length,2);assert(entries.every(r=>r.annotationId===region.annotation.id&&r.modelProteinId===region.modelProteinId&&r.name&&r.fullName));
 for(const edit of [r=>r.id='unknown',r=>r.modelProteinId='source-b',r=>r.correspondence[0].canonical++]){const bad=structuredClone(region);edit(bad);assert.throws(()=>siteResidueEntries(data,ledger,loc,bad));}
 const other=endpointLocation(data,ledger,'pair-1','source-b');assert.equal(other.regions.length,0);assert.deepEqual(siteResidueEntries(data,ledger,other,null),[]);assert.equal(JSON.stringify(data),before);
});

test('large constructs retain every selectable residue while projected labels are bounded and selection is retained',()=>{
 const rows=Array.from({length:708},(_,index)=>({index,name:'ALA',letter:'A',fullName:'Alanine',chain:'B',author:index+17,authorLabel:'B:'+(index+17),paper:index+17,canonical:index+41})),before=JSON.stringify(rows),indices=[];
 for(let page=0;page<30;page++)indices.push(...siteResiduePage(rows,{page}).items.map(r=>r.index));assert.equal(indices.length,708);assert.equal(new Set(indices).size,708);
 const labels=siteLabelCandidates(rows,707);assert(labels.length<=12);assert.equal(labels[0].index,707);assert.equal(new Set(labels.map(r=>r.index)).size,labels.length);
 assert.equal(siteResiduePage(rows,{query:'B:724'}).items[0].index,707);assert.equal(siteResiduePage(rows,{query:'canonical 748'}).items[0].index,707);assert.equal(siteResiduePage(rows,{query:'unknown-homolog'}).total,0);assert.equal(siteLabelCandidates(rows.slice(0,2)).length,2);assert.equal(JSON.stringify(rows),before);
});

test('residue evidence text cannot promote a construct, secondary cut, projection or prediction into an atomic contact',()=>{
 const c=residueEvidenceClaim({category:'assayed_construct'});assert.match(c.title,/membership/);assert.match(c.claim,/does not identify/);
 assert.match(residueEvidenceClaim({category:'measured_reactive_bond'}).claim,/not a solved/);
 assert.match(residueEvidenceClaim({category:'measured_secondary_cleavage'}).claim,/not assigned as an inhibitory/);
 assert.match(residueEvidenceClaim({category:'functional_region',annotation:{status:'conditional_reference_projection'}}).title,/Conditional/);
 assert.match(residueEvidenceClaim({category:'prediction'}).claim,/not measured/);
});
