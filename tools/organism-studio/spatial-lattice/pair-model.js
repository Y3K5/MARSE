// Side-by-side display of literal source endpoints; never a docking pose or a new relationship.
import {supportedAcrossPairs,foldRadius,displayLimits} from './exploration-model.js';
export const pairLimits={min:1,max:4,default:2.4};
export function validatePairFocus(data,state){
 if(!state||typeof state!=='object'||Array.isArray(state)||Object.keys(state).some(k=>!['record','magnification'].includes(k))||!Number.isFinite(state.magnification)||state.magnification<pairLimits.min||state.magnification>pairLimits.max)throw Error('Invalid large-pair display settings.');
 const record=supportedAcrossPairs(data).find(r=>r.id===state.record);if(!record)throw Error('No supported literal cross-microbe protein pair for this display.');return record;
}
export function pairSceneData(data,state){
 const record=validatePairFocus(data,state),ids=new Set([record.a,record.b]);
 return {...data,proteins:data.proteins.filter(p=>ids.has(p.id)),contexts:[],relationships:[record],display_field:'source-pair'};
}
export function pairLayout(data,state,foldSize,spacing){
 const record=validatePairFocus(data,state);
 if(!Number.isFinite(foldSize)||foldSize<displayLimits.foldMin||foldSize>displayLimits.foldMax||!Number.isFinite(spacing)||spacing<displayLimits.spacingMin||spacing>displayLimits.spacingMax)throw Error('Invalid pair display scale.');
 const radius=foldRadius(foldSize)*state.magnification,separation=2*radius+14*spacing+8;
 return Object.fromEntries([record.a,record.b].map((id,i)=>{const p=data.proteins.find(v=>v.id===id),position=[(i-.5)*separation,0,0];return [id,{id,name:p.name,species:p.species,kind:'protein',hasShape:!!data.shape_traces[id],position,groupCenter:[...position],displayRadius:radius,clearancePitch:separation}];}));
}
