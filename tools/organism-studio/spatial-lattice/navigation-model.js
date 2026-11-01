// Viewing/navigation only. Never changes target membership, source evidence or relationships.
import {resultPage,sceneData} from './subtractive-model.js';
import {supportedAcrossPairs} from './exploration-model.js';
export function proteinNavigationItem(data,id){
 const p=data.proteins.find(p=>p.id===id);if(!p)return null;
 const r=p.subtractive_record,sp=data.species.find(s=>s.id===p.species),hasShape=!!data.shape_traces[id];
 return {id,kind:'protein',name:p.name,scope:r?'Retained result':'Source endpoint',identity:[r?.strain||sp?.short||sp?.name||p.species,p.accession||r?.key||'Reference unresolved'].join(' · '),hasShape,status:hasShape?(p.retained_evidence?.structure?.label||'Retained source fold'):'Coordinate hold · no fold assigned'};
}
export function navigationItems(data,options={},pair=false){
 if(pair)return supportedAcrossPairs(data).map(r=>({id:r.id,kind:'pair',name:[r.a,r.b].map(id=>data.proteins.find(p=>p.id===id).name).join(' ↔ '),scope:'Published-source relationship',identity:[r.a,r.b].map(id=>proteinNavigationItem(data,id).identity).join(' ↔ '),status:r.level+' · source strain/material limits apply',hasShape:[r.a,r.b].some(id=>data.shape_traces[id])}));
 return (data.subtractive?resultPage(data,options).items:data.proteins).map(p=>proteinNavigationItem(data,p.id));
}
export function navigationCursor(items,id){const index=items.findIndex(v=>v.id===id);return {index,total:items.length,previous:index>0?items[index-1].id:null,next:index<0?items[0]?.id||null:items[index+1]?.id||null};}
export function validateNavigationProtein(data,options,id){
 if(typeof id!=='string'||!(data.subtractive?sceneData(data,options).proteins:data.proteins).some(p=>p.id===id))throw Error('Saved protein selection is outside the visible reference/filter scope.');return id;
}
