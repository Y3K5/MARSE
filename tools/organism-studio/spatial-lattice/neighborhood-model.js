// Local evidence exploration only. No inferred partner, site mapping or target adoption.
import {isSupportedProteinPair,sourceReferences,reviewEntry} from './review-model.js';
import {displayLimits,foldRadius} from './exploration-model.js';
import {endpointLocation} from './interaction-model.js';
export function validateNeighborhood(data,state){
 if(!state||typeof state!=='object'||Array.isArray(state)||Object.keys(state).some(k=>k!=='protein')||typeof state.protein!=='string'||!data.proteins.some(p=>p.id===state.protein))throw Error('Unknown protein neighborhood; no endpoint is substituted.');return state;
}
export function proteinNeighborhood(data,state){
 validateNeighborhood(data,state);const seed=data.proteins.find(p=>p.id===state.protein),references=sourceReferences(data,seed.id),sourceIds=seed.subtractive_record?references.map(l=>l.b):[seed.id],sources=new Set(sourceIds);
 const incident=data.relationships.filter(r=>sources.has(r.a)||sources.has(r.b)),records=incident.filter(r=>isSupportedProteinPair(data,r)),ids=new Set([seed.id,...sourceIds,...records.flatMap(r=>[r.a,r.b])]);
 const proteins=[seed,...data.proteins.filter(p=>p.id!==seed.id&&ids.has(p.id))],holds=reviewEntry(data,seed.id)?.holds||[];
 return {seed,sourceIds,references,records,otherRecords:incident.filter(r=>!records.includes(r)&&r.level!=='candidate / unresolved identity'),proteins,holds,partners:[...new Set(records.flatMap(r=>[r.a,r.b]).filter(id=>!sources.has(id)))],scope:seed.subtractive_record?'Retained reference with separate source context':'Literal published-source endpoint'};
}
export function neighborhoodScene(data,state){const n=proteinNeighborhood(data,state);return {...data,proteins:n.proteins,contexts:[],relationships:n.records,display_field:'protein-neighborhood'};}
export function neighborhoodLayout(data,state,size,spacing,mode='cluster'){
 if(!Number.isFinite(size)||size<displayLimits.foldMin||size>displayLimits.foldMax||!Number.isFinite(spacing)||spacing<displayLimits.spacingMin||spacing>displayLimits.spacingMax||!['cluster','layers','orbit'].includes(mode))throw Error('Invalid neighborhood display settings.');
 const n=proteinNeighborhood(data,state),radius=foldRadius(size)*1.65,pitch=2*radius+12*spacing+8,count=n.proteins.length-1,ring=count<2?pitch:Math.max(pitch,pitch/(2*Math.sin(Math.PI/count))),out={};
 n.proteins.forEach((p,i)=>{const angle=(i-1)/Math.max(1,count)*Math.PI*2+({cluster:.7,layers:.2,orbit:.4}[mode]),position=i?[Math.cos(angle)*ring,((i-1)%3-1)*pitch*.3,Math.sin(angle)*ring]:[0,0,0];out[p.id]={id:p.id,name:p.name,species:p.species,kind:'protein',hasShape:!!data.shape_traces[p.id],position,groupCenter:[0,0,0],displayRadius:radius,clearancePitch:pitch};});return out;
}
export function neighborhoodSites(data,ledger,state){const n=proteinNeighborhood(data,state),sites=[];for(const record of n.records)for(const endpoint of new Set([record.a,record.b])){const location=endpointLocation(data,ledger,record.id,endpoint);for(const region of location.regions)sites.push({record:record.id,endpoint,region,location});}return sites;}
export function validateNeighborhoodSelection(data,state,id){if(!neighborhoodScene(data,state).proteins.some(p=>p.id===id))throw Error('Saved protein is outside this qualified neighborhood.');return id;}

// Source-coordinate cameras use a different frame and must not be replayed into a display field.
export function neighborhoodReturnView(preset,context,camera){
 return preset==='microbes'&&context.kind==='network'?{position:[...camera.position],target:[...camera.target]}:null;
}
