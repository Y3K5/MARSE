import assert from 'node:assert/strict';
import test from 'node:test';
import {createHash} from 'node:crypto';
import {AtomicStore,parseAtomicPdb,validateAtomicPayload,atomicNeighbors,atomicEnvironment} from './atomic-model.js';
const hash=b=>createHash('sha256').update(b).digest('hex');
function line(serial,name,residue,x,element,ins=''){
 const chars=Array(80).fill(' '),put=(start,width,value)=>chars.splice(start,width,...String(value).padStart(width));
 put(0,6,'ATOM  ');put(6,5,serial);put(12,4,name);put(17,3,'ALA');put(21,1,'A');put(22,4,residue);put(26,1,ins);put(30,8,x.toFixed(3));put(38,8,'0.000');put(46,8,'0.000');put(54,6,'1.00');put(60,6,'20.00');put(76,2,element);return chars.join('');
}
const pdb=[line(1,'CA',8,0,'C','A'),line(2,'C',8,1.4,'C','A'),line(3,'N',9,2.7,'N'),line(4,'CA',9,4,'C')].join('\n'),pin=hash(pdb),data={shape_traces:{invented:{source_view_sha256:pin,chains:[{chain:'A',ca:[[8,0,0,0],[9,4,0,0]]}]}}},ledger={entries:{invented:{source_view_sha256:pin,chains:[{chain:'A',residues:[[8,'ALA','A'],[9,'ALA','']]}]}}},payload={schema:'marse-atomic-context-v1',protein_id:'invented',source_view_sha256:pin,scientific_execution:false,context:{source_url:'https://example.org/structure',kind:'experimental'},atoms:parseAtomicPdb(pdb),bonds:[{a:1,b:2,order:'sing',kind:'polymer_link',distance:1.3}]};
const valid=p=>validateAtomicPayload(p,pdb,data,ledger,'invented');
test('exact retained atoms preserve insertion/chain/name/element/coordinates and topology is separate',()=>{
 const before=JSON.stringify({payload,data,ledger});valid(payload);assert.equal(atomicNeighbors(payload,1)[0].atom.name,'N');assert.deepEqual(atomicEnvironment(payload,1,1.4),[0,1,2]);assert.equal(atomicNeighbors(payload,0).length,0);assert.equal(JSON.stringify({payload,data,ledger}),before);
 for(const mutate of [p=>p.atoms[1].name='O',p=>p.atoms[1].xyz[0]+=.01,p=>p.atoms[1].element='N',p=>p.atoms[0].insertion='',p=>p.atoms[0].chain='B',p=>p.atoms.pop(),p=>p.source_view_sha256='0'.repeat(64),p=>p.protein_id='homolog',p=>p.context.kind='measured',p=>p.context.source_url='javascript:alert(1)',p=>p.context.source_url='https://private@example.org',p=>p.scientific_execution=true]){const p=structuredClone(payload);mutate(p);assert.throws(()=>valid(p));}
});
test('unsupported, duplicate, self and geometry-altered links are held',()=>{
 for(const mutate of [p=>p.bonds[0].a=-1,p=>p.bonds[0].b=4,p=>p.bonds[0].b=1,p=>p.bonds[0].kind='nearby_contact',p=>p.bonds[0].order='guessed',p=>p.bonds[0].distance+=.01,p=>p.bonds.push(p.bonds[0])]){const p=structuredClone(payload);mutate(p);assert.throws(()=>valid(p));}
 assert.throws(()=>atomicEnvironment(payload,99));assert.throws(()=>atomicEnvironment(payload,1,30));
});
function fixture(edit=()=>{}){
 const p=JSON.stringify(payload),name=hash(p)+'.json',pdbName=pin+'.pdb',index={schema:'marse-atomic-index-v1',scientific_execution:false,source_snapshot_sha256:'snapshot',residue_ledger_sha256:'residues',entries:{invented:{file:name,sha256:hash(p),pdb_file:pdbName,source_view_sha256:pin}}};edit(index);const raw=JSON.stringify(index),receipt=JSON.stringify({schema:'marse-atomic-import-v1',index_sha256:hash(raw),files:{[name]:hash(p),[pdbName]:pin}}),files={'index.json':raw,'IMPORT.json':receipt,[name]:p,[pdbName]:pdb};let calls=0;return {files,store:new AtomicStore(data,ledger,'snapshot',async url=>{calls++;return new Response(files[url.replace('runtime/atomic/','')]??'',{status:files[url.replace('runtime/atomic/','')]?200:404});},'residues'),calls:()=>calls};
}
test('reader validates byte pins, snapshot, ledger, safe paths, exact fold and cache; no substitution',async()=>{
 const f=fixture();await f.store.load('invented');const count=f.calls();await f.store.load('invented');assert.equal(f.calls(),count);await assert.rejects(f.store.load('homolog'));
 for(const mutate of [i=>i.source_snapshot_sha256='old',i=>i.residue_ledger_sha256='old',i=>i.entries.invented.file='../private.json',i=>i.entries.invented.pdb_file='../model.pdb',i=>i.entries.invented.source_view_sha256='bad'])await assert.rejects(fixture(mutate).store.load('invented'));
 const bad=fixture();bad.files[pin+'.pdb']+='\n';await assert.rejects(bad.store.load('invented'));
 const badIndex=fixture();badIndex.files['index.json']+=' ';await assert.rejects(badIndex.store.load('invented'));
});

test('incomplete PDB numeric fields and nonfinite metadata are held instead of becoming zero coordinates',()=>{
 const first=pdb.split('\n')[0];
 for(const [offset,width,value] of [[6,5,''],[22,4,''],[30,8,''],[38,8,''],[46,8,''],[54,6,'NaN'],[60,6,'Inf'],[12,4,''],[17,3,'']]){
  const bad=first.slice(0,offset)+value.padStart(width)+first.slice(offset+width);assert.throws(()=>parseAtomicPdb(bad));
 }
});

test('private atomic successor activation is pinned and fails closed; legacy import stays readable',async()=>{
 const f=fixture(),indexHash=hash(f.files['index.json']),pointer={schema:'marse-atomic-active-v1',directory:'atomic-coverage-v1',index_sha256:indexHash,source_snapshot_sha256:'snapshot',residue_ledger_sha256:'residues'};
 function store(p=pointer){return new AtomicStore(data,ledger,'snapshot',async url=>{
  if(url==='runtime/ATOMIC_ACTIVE.json')return new Response(JSON.stringify(p));
  const name=url.replace('runtime/atomic-coverage-v1/','');return new Response(f.files[name]??'',{status:f.files[name]?200:404});
 },'residues');}
 assert.equal((await store().load('invented')).atoms.length,4);
 for(const mutate of [p=>p.directory='../private',p=>p.directory='https://example.org',p=>p.directory='atomic/x',p=>p.index_sha256='0'.repeat(64),p=>p.source_snapshot_sha256='stale',p=>p.residue_ledger_sha256='stale',p=>p.schema='untrusted']){
  const p=structuredClone(pointer);mutate(p);await assert.rejects(store(p).load('invented'));
 }
 await f.store.load('invented');assert.equal(f.store.root,'atomic');
});

test('explicit coordinate holds do not substitute a matching-looking model',async()=>{
 const f=fixture(index=>{index.holds=[{protein_id:'missing',reason:'No retained matching full-coordinate structure'}];});
 await assert.rejects(f.store.load('missing'),/No retained matching full-coordinate structure.*No related protein is substituted/);
});
