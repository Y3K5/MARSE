import {isSupportedProteinPair} from './review-model.js';
// Display organization only. Protein identities, coordinates and relationships stay in the packet.
export const displayLimits={foldMin:1,foldMax:4.5,foldDefault:3.2,spacingMin:.7,spacingMax:3,spacingDefault:1.5};
export const foldRadius=size=>4.6*size;
export const haloRadius=size=>5.2*size;

export function retainedGroups(data){
 if(!data.subtractive)return [];
 const review=new Set(data.subtractive.shortlist_keys);
 return data.subtractive.species.map(s=>{
  const rows=data.proteins.filter(p=>p.subtractive_record?.species===s.id);
  const forms=rows.filter(p=>p.retained_evidence?.structure?.status!=='missing');
  const folds=forms.filter(p=>!!data.shape_traces[p.id]).length;
  return {...s,nativeSpecies:rows[0]?.species||null,forms:forms.length,folds,holds:forms.length-folds,review:rows.filter(p=>review.has(p.subtractive_record.key)).length};
 });
}

export function explorationChoices(data){
 const groups=retainedGroups(data),review=groups.filter(g=>g.review),choices=[],reviewKeys=new Set(data.subtractive?.shortlist_keys||[]),reviewGaps=data.proteins.filter(p=>reviewKeys.has(p.subtractive_record?.key)&&p.retained_evidence?.structure?.status==='missing').length;
 if(review.length)choices.push({id:'review',label:review.length===1?(review[0].strain||review[0].label||review[0].id)+' review':'Review queue',folds:data.proteins.filter(p=>reviewKeys.has(p.subtractive_record?.key)&&data.shape_traces[p.id]).length,records:review.reduce((n,g)=>n+g.review,0)});
 for(const g of groups.filter(g=>!g.review))choices.push({id:'forms:'+g.id,label:g.label,folds:g.folds,holds:g.holds,records:g.forms});
 choices.unshift({id:'together',label:'All microbes + evidence',folds:groups.reduce((n,g)=>n+g.folds,0),holds:groups.reduce((n,g)=>n+g.holds,0)+reviewGaps,sourceFolds:sourceNetworkData(data).proteins.filter(p=>data.shape_traces[p.id]).length});
 choices.push({id:'all',label:'All retained references',folds:groups.reduce((n,g)=>n+g.folds,0),holds:groups.reduce((n,g)=>n+g.holds,0),records:groups.reduce((n,g)=>n+g.forms,0)});
 return choices;
}

export function explorationOptions(data,id){
 if(id==='together')return {cohort:'forms',organism:'all',sourceContext:true};
 if(id==='review')return {cohort:'shortlist',organism:'all'};
 if(id==='all')return {cohort:'forms',organism:'all'};
 if(id.startsWith('forms:')&&retainedGroups(data).some(g=>g.id===id.slice(6)))return {cohort:'forms',organism:id.slice(6)};
 throw Error('Unknown retained result view.');
}

export function sourceNetworkData(data){
 if(!data.subtractive)return data;
 const ids=new Set(data.subtractive.source_context_ids),proteins=data.proteins.filter(p=>ids.has(p.id)),contexts=data.contexts;
 for(const p of contexts)ids.add(p.id);
 return {...data,proteins,contexts,relationships:data.relationships.filter(r=>ids.has(r.a)&&ids.has(r.b))};
}

// Keep supported source pairs literal: cell/class, phenotype, preparation-only and candidate records are excluded.
export function supportedAcrossPairs(data){
 const by=new Map(data.proteins.map(p=>[p.id,p]));
 return data.relationships.filter(r=>isSupportedProteinPair(data,r)&&by.get(r.a).species!==by.get(r.b).species);
}

export function fieldLayout(data,spacing,mode,size,centers){
 const pitch=2*haloRadius(size)+12*spacing,groups=[];
 for(const species of Object.keys(centers)){
  const members=[...data.proteins,...data.contexts].filter(p=>p.species===species).sort((a,b)=>Number(!!b.subtractive_record)-Number(!!a.subtractive_record)||String(a.subtractive_record?.primary_group||'source').localeCompare(String(b.subtractive_record?.primary_group||'source'))||String(a.id).localeCompare(String(b.id)));
  if(!members.length)continue;
  const n=members.length,cols=mode==='cluster'?Math.ceil(Math.cbrt(n)):Math.ceil(Math.sqrt(n*(data.display_field==='combined'?.8:2))),depth=mode==='cluster'?Math.ceil(n/(cols*cols)):Math.ceil(n/cols),shelfCount=Math.min(3,n),perShelf=Math.ceil(n/shelfCount);
  const ringRadius=perShelf<2?pitch:Math.max(pitch,pitch/(2*Math.sin(Math.PI/perShelf)));
  const offsets=members.map((p,i)=>{
   if(mode==='cluster')return [(i%cols-(cols-1)/2)*pitch,(Math.floor(i/(cols*cols))-(depth-1)/2)*pitch,(Math.floor(i/cols)%cols-(cols-1)/2)*pitch];
   if(mode==='layers')return [(i%cols-(cols-1)/2)*pitch,((Math.floor(i/cols)%3)-1)*pitch*.42,(Math.floor(i/cols)-(depth-1)/2)*pitch];
   const shelf=Math.floor(i/perShelf),angle=(i%perShelf)/perShelf*Math.PI*2+shelf*.22;
   return [Math.cos(angle)*ringRadius,(shelf-(shelfCount-1)/2)*pitch*1.35,Math.sin(angle)*ringRadius];
  });
  groups.push({species,members,offsets,base:centers[species]});
 }
 if(!groups.length)return {};
 const extent=Math.max(...groups.flatMap(g=>g.offsets.map(o=>Math.hypot(...o))));
 let minDistance=Infinity;
 for(let i=0;i<groups.length;i++)for(let j=0;j<i;j++)minDistance=Math.min(minDistance,Math.hypot(...groups[i].base.map((v,k)=>v-groups[j].base[k])));
 const scale=groups.length>1?(2*extent+3*pitch)/minDistance:0;
 // Combined display packs species fields along X using their actual extents; no spherical empty moat.
 const packed=new Map();
 if(data.display_field==='combined'){
  let cursor=0;
  for(const g of groups){const lo=Math.min(...g.offsets.map(o=>o[0])),hi=Math.max(...g.offsets.map(o=>o[0]));packed.set(g.species,[cursor-lo,0,0]);cursor+=hi-lo+3*pitch;}
  const middle=(cursor-3*pitch)/2;for(const v of packed.values())v[0]-=middle;
 }
 const mean=[0,1,2].map(k=>groups.reduce((n,g)=>n+g.base[k],0)/groups.length),out={};
 for(const g of groups){
  const center=packed.get(g.species)||g.base.map((v,k)=>(v-mean[k])*scale);
  g.members.forEach((p,i)=>{out[p.id]={id:p.id,name:p.name,species:g.species,kind:data.contexts.some(c=>c.id===p.id)?(['class','protein_class'].includes(p.kind)?'class':['cell','whole_cell'].includes(p.kind)?'cell':'context'):'protein',hasShape:!!data.shape_traces[p.id],position:center.map((v,k)=>v+g.offsets[i][k]),groupCenter:[...center],displayRadius:foldRadius(size),clearancePitch:pitch};});
 }
 return out;
}

// Preserve unassigned source/context markers without inventing a fourth organism.
export function overviewPointGroups(layout){
 const groups=new Map();for(const node of Object.values(layout)){if(!groups.has(node.species))groups.set(node.species,[]);groups.get(node.species).push(node);}
 return [...groups].map(([species,nodes])=>({species,nodes}));
}
