import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {MORPHOLOGY,morphologySVG,morphologyPanel} from '../dist/pocket-core/pocket-morphology.mjs';
import {TAXA,observation} from '../dist/pocket-core/pocket-science.mjs';
import {SEQUENCES} from '../dist/pocket-core/pocket-sequences.mjs';

test('every saved taxon has a separate qualitative form and source context',()=>{
 assert.deepEqual(Object.keys(MORPHOLOGY),TAXA.map(t=>t.id));
 for(const t of TAXA){
  const m=MORPHOLOGY[t.id];assert.match(m.source,/^https:\/\//);assert.ok(m.context.length>40);
  assert.match(morphologySVG(t.id),new RegExp(t.name.replace(/[.*+?^${}()|[\]\\]/g,'\\$&')));
  assert.ok(morphologyPanel(t.id).includes(m.source));
  assert.match(morphologyPanel(t.id),/no physical scale/);
 }
 assert.equal(MORPHOLOGY.pg.form,'short-rod');assert.equal(MORPHOLOGY.td.form,'spirochete');
 assert.match(MORPHOLOGY.ac.observation,/not a species identification/);
 assert.equal((morphologyPanel('all').match(/data-morphology=/g)||[]).length,9);
 assert.throws(()=>morphologySVG('unknown'),/Unknown morphology identity/);
});
test('morphology illustrations do not assign new spatial evidence or mutate saved fields',async()=>{
 const path=new URL('../src/pocket-input.json',import.meta.url),before=await readFile(path);
 for(const t of TAXA){morphologyPanel(t.id);assert.equal(observation(t.id,'healthy').resolved,false);}
 assert.equal(observation('td','periodontitis').resolved,false);
 assert.deepEqual(await readFile(path),before);
});

test('storyboard keeps bone loss separate from a calibrated disease clock',()=>{
 assert.equal(SEQUENCES.response[5].source,'boneDefects');
 assert.match(SEQUENCES.response[4].text,/not an inevitable outcome/);
 assert.equal(SEQUENCES.response[5].status,'Illustrative');
});
