'use strict';
// Optional reference checks; no dependency on the MARSE Python solver.
const assert=require('node:assert/strict');
const fs=require('node:fs'),path=require('node:path');
const K=require('./reference.cjs');
const read=name=>fs.readFileSync(path.join(__dirname,name));
const protocol=JSON.parse(read('protocol.json'));
const tolerance=protocol.tolerances;
function inspect(run){
  assert.ok(run.audit.finite&&run.audit.minimum>=0);
  assert.ok(Math.abs(run.audit.biomassResidual)<tolerance.ledger_absolute);
  run.audit.resourceResiduals.forEach(v=>assert.ok(Math.abs(v)<tolerance.ledger_absolute));
  assert.ok(run.ledger.maximumResourceOverdraw<=tolerance.resource_overdraw_absolute);
}
const low=K.simulate({duration:.01,resourceInitial:[.0001,0],taxa:[
  {id:'synthetic-a',rate:1e6,mortality:0,initial:.01,costs:[.5,0]},
  {id:'synthetic-b',rate:1e6,mortality:0,initial:.01,costs:[.5,0]}
]});inspect(low);
assert.ok(Math.abs(low.ledger.birth*.5-low.ledger.consumed[0])<tolerance.ledger_absolute);
const empty=K.simulate({resourceInitial:[0,0]});inspect(empty);assert.equal(empty.ledger.birth,0);
const dying=K.simulate({taxa:[{id:'synthetic-a',rate:0,mortality:200,initial:.1,costs:[1,0]}]});inspect(dying);
const boundary=K.simulate({boundaryRate:1e5,boundaryReservoir:[0,0],taxa:[{id:'synthetic-a',rate:500,mortality:200,initial:.1,costs:[.5,.5]}]});inspect(boundary);
const forward=K.simulate(),reverse=K.simulate({taxa:[...forward.config.taxa].reverse()});inspect(forward);inspect(reverse);assert.deepEqual(forward.state,reverse.state);
const before=K.initialize(K.configuration()),inserted=K.initialize(K.configuration({taxa:[...forward.config.taxa,{id:'synthetic-unrelated',rate:0,mortality:0,initial:.01,costs:[1,0]}]}));
for(const id of Object.keys(before.biomass))assert.deepEqual(before.biomass[id],inserted.biomass[id]);
assert.deepEqual(forward.state,K.simulate().state);assert.throws(()=>K.simulate({dt:1}),/positivity bound/);
const c=K.configuration({nx:9,ny:17}),a=Float64Array.from({length:c.nx*c.ny},(_,i)=>i%7*.01);
const fluxResidual=K.integral(K.diffuse(a,c,.01,.001),c)-K.integral(a,c);assert.ok(Math.abs(fluxResidual)<tolerance.ledger_absolute);
const diff=protocol.diffusion_convergence,decay=Math.exp(-2*Math.PI**2*diff.diffusion*diff.duration);
const diffusion=diff.grids.map(([nx,ny])=>{
  const run=K.simulate({nx,ny,dt:diff.dt,duration:diff.duration,diffusion:diff.diffusion,smoothInitial:true,taxa:[{id:'synthetic-a',rate:0,mortality:0,initial:1,costs:[1,0]}]});inspect(run);
  let errors=0,normalizer=0;run.state.biomass['synthetic-a'].forEach((v,i)=>{
    const exact=1+.2*Math.cos(Math.PI*(i%nx+.5)/nx)*Math.cos(Math.PI*(Math.floor(i/nx)+.5)/ny)*decay;
    errors+=(v-exact)**2;normalizer+=exact**2;
  });return {nx,ny,dt:diff.dt,error:Math.sqrt(errors/normalizer)};
});
const t=protocol.time_convergence,reference=K.simulate({nx:t.grid[0],ny:t.grid[1],dt:t.reference_dt,duration:t.duration,smoothInitial:true});inspect(reference);
const time=t.steps.map(dt=>{const r=K.simulate({nx:t.grid[0],ny:t.grid[1],dt,duration:t.duration,smoothInitial:true});inspect(r);return {dt,error:Math.abs(r.final.biomass-reference.final.biomass)/reference.final.biomass};});
const g=protocol.reaction_transport_grid,ref=K.simulate({nx:g.reference_grid[0],ny:g.reference_grid[1],dt:g.dt,duration:g.duration,smoothInitial:true});inspect(ref);
const grid=g.grids.map(([nx,ny])=>{const r=K.simulate({nx,ny,dt:g.dt,duration:g.duration,smoothInitial:true});inspect(r);return {nx,ny,dt:g.dt,error:Math.abs(r.final.biomass-ref.final.biomass)/ref.final.biomass};});
for(const rows of [diffusion,time,grid])for(let i=1;i<rows.length;i++)assert.ok(rows[i].error<rows[i-1].error);
assert.ok(diffusion.at(-1).error<.003&&time.at(-1).error<.005&&grid.at(-1).error<.005);
console.log(JSON.stringify({engine:K.ENGINE,sourceHash:K.digest(read('reference.cjs')),protocolHash:K.digest(read('protocol.json')),
  statuses:{implementation:'implemented',execution:'executed_synthetic',gate:'passed',evidence:'synthetic_only',approval:'review_required'},
  checks:{resourceSupportedGrowth:true,zeroResourceBirth:empty.ledger.birth,nonnegativeStressStates:true,identitySeedInvariant:true,permutationInvariant:true,reproducible:true,fluxResidual},
  convergence:{diffusion,time,grid},calibration:'not_run',biologicalEvaluation:'not_run',marseEngineValidation:'not_established_by_this_reference'},null,2));
