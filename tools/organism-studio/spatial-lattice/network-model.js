import {resultLayout} from './subtractive-model.js';
import {displayLimits} from './exploration-model.js';
export const networkCenters={};
export const neighborhoodViews={cluster:{name:'Clustered neighborhood',description:'Retained folds in spatial groups',centers:networkCenters,direction:[.18,.45,1.6],radius:110,focus:[0,0,0]},layers:{name:'Layered neighborhood',description:'Staggered source groups',centers:{},direction:[.2,1.8,1.2],radius:110,focus:[0,0,0]},orbit:{name:'Orbital neighborhood',description:'Source groups on circular shelves',centers:{},direction:[.22,.68,1.6],radius:110,focus:[0,0,0]}};
export function configureNetworkSpecies(data){const ids=[...data.species.map(s=>s.id)];if([...data.proteins,...data.contexts].some(p=>p.species==='unknown')&&!ids.includes('unknown'))ids.push('unknown');ids.forEach((id,i)=>{const a=i/Math.max(1,ids.length)*Math.PI*2;networkCenters[id]=[Math.cos(a)*45,(i%3-1)*14,Math.sin(a)*45];neighborhoodViews.layers.centers[id]=[(i%2?1:-1)*30,(i-(ids.length-1)/2)*34,Math.cos(a)*20];neighborhoodViews.orbit.centers[id]=[Math.cos(a)*52,(i%3-1)*13,Math.sin(a)*52];});}
export function neighborhoodView(id='cluster'){const view=neighborhoodViews[id];if(!Object.hasOwn(neighborhoodViews,id)||!view)throw Error('Unknown microbial neighborhood layout.');return view;}
export const networkStyles={
 binding:{label:'Binding / protein interaction',color:'#7bd8cf',pattern:'solid'},
 inhibition:{label:'Biochemical inhibition',color:'#ee9d96',pattern:'solid'},
 association:{label:'Association / assembly',color:'#b9abe5',pattern:'solid'},
 processing:{label:'Processing / domain relationship',color:'#d8ba80',pattern:'dash'},
 context:{label:'Cell, phenotype or assembly context',color:'#adc0cd',pattern:'dash'},
 reported:{label:'Context-dependent binding preparation',color:'#93a9e4',pattern:'dash'},
 candidate:{label:'Candidate / unresolved identity',color:'#d6ab64',pattern:'dot'}
};
const levels={'direct assay':'binding','biochemical interaction':'association','direct biochemical':'inhibition','biochemical complex':'association','assembly / microscopy':'association','structural and biochemical':'association','experimental structural assembly':'association','biochemical processing':'processing','domain / processing':'processing','assembly context':'context','indirect / phenotype':'context','direct cell-binding / receptor unresolved':'context','community phenotype':'context','binding preparation / context dependent':'reported','candidate / unresolved identity':'candidate'};
export function recordStyle(record){const key=levels[record.level];if(!key)throw Error('Unrecognized source evidence class: '+record.level);return {key,...networkStyles[key]};}
export function networkRecords(data,{scope='all',candidates=false}={}){
 if(!['all','across','direct'].includes(scope)||typeof candidates!=='boolean')throw Error('Invalid network filter.');
 const by=new Map([...data.proteins,...data.contexts].map(p=>[p.id,p]));
 return data.relationships.filter(r=>{const a=by.get(r.a),b=by.get(r.b),style=recordStyle(r);if(!a||!b)throw Error('Missing literal source endpoint.');if(!candidates&&style.key==='candidate')return false;if(scope==='across'&&a.species===b.species)return false;if(scope==='direct'&&!['direct assay','direct biochemical'].includes(r.level))return false;return true;});
}
export function networkLayout(data,spacing=1,mode='cluster',foldSize=1){
 const centers=neighborhoodView(mode).centers;if(!Number.isFinite(foldSize)||foldSize<displayLimits.foldMin||foldSize>displayLimits.foldMax)throw Error('Invalid fold display size.');
 if(!Number.isFinite(spacing)||spacing<displayLimits.spacingMin||spacing>displayLimits.spacingMax)throw Error('Invalid network spacing.');const result={};
 return resultLayout(data,spacing,mode,foldSize,centers);
}

export function edgePath(record,layout,index=0){
 const a=layout[record.a]?.position,b=layout[record.b]?.position;if(!a||!b)throw Error('Missing network endpoint position.');
 if(record.a===record.b)return [a,[a[0]+7,a[1]+8,a[2]+3],[a[0]-6,a[1]+9,a[2]-6],a];
 const arc=(index%5-2)*2.1,mid=a.map((v,j)=>(v+b[j])/2+(j===1?arc:j===2?3+(index%3)*2:0));return [a,mid,b];
}
export function networkStats(data,records){
 const all=[...data.proteins,...data.contexts],by=new Map(all.map(p=>[p.id,p])),adj=new Map(),types={};let across=0,contexts=0;
 for(const r of records){if(!by.has(r.a)||!by.has(r.b))throw Error('Missing literal endpoint.');for(const [a,b]of [[r.a,r.b],[r.b,r.a]]){if(!adj.has(a))adj.set(a,new Set());adj.get(a).add(b);}const key=recordStyle(r).key;types[key]=(types[key]||0)+1;if(by.get(r.a).species!==by.get(r.b).species)across++;if(data.contexts.some(p=>p.id===r.a||p.id===r.b))contexts++;}
 const seen=new Set(),components=[];for(const id of adj.keys()){if(seen.has(id))continue;const todo=[id],component=[];seen.add(id);while(todo.length){const n=todo.pop();component.push(n);for(const p of adj.get(n))if(!seen.has(p)){seen.add(p);todo.push(p);}}components.push(component.sort());}
 return {records:records.length,across,within:records.length-across,types,linkedNodes:[...adj.keys()],isolatedNodes:all.filter(p=>!adj.has(p.id)).map(p=>p.id),components,contextRecords:contexts};
}
