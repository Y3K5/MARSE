import type {PocketId} from './pocket-geometry';
import {PRESETS,type Health,type Scale} from './pocket-science';

const primary:PocketId[]=['enamel','dentin','pulp','gingiva','supragingival','cementum','plaque','lumen','epithelium','pdl','bone','cheek','tongue'];
export const sceneLabelIds:PocketId[]=[...primary,'connective','vessels','nerve'];
// Schematic search targets; labels snap to actual exported vessel/nerve vertices.
export const additionalLabelTargets:Partial<Record<PocketId,[number,number,number]>>={vessels:[3.1,-8,0],nerve:[3.25,-8,0]};
const inside:PocketId[]=['dentin','pulp','cementum','plaque','lumen','epithelium','pdl','connective','vessels','nerve'];
export const shortLabels:Partial<Record<PocketId,string>>={pulp:'Pulp',supragingival:'Surface plaque',plaque:'Pocket plaque',lumen:'Fluid · widened',epithelium:'Lining',pdl:'Ligament',bone:'Bone',connective:'Connective',vessels:'Vessels',nerve:'Nerve',cheek:'Cheek',tongue:'Tongue'};
export function visibleLabels(width:number,cutaway:boolean,scale:Scale,selected:PocketId){
 const available=sceneLabelIds.filter(id=>!(!cutaway&&inside.includes(id))&&!(['cheek','tongue'].includes(id)&&scale!=='tooth'));
 if(width<560){
  const essential=['gingiva','plaque','bone',available.includes(selected)?selected:'enamel'] as PocketId[];
  return [...new Set(essential)].filter(id=>available.includes(id));
 }
 // A few anchor labels, with every tissue reachable through the inspector.
 return available.filter(id=>['enamel','gingiva','plaque','epithelium','pdl','bone',selected].includes(id));
}
export function sceneCopy(scale:Scale,health:Health,cutaway:boolean){
 const form=cutaway?'quarter cutaway':'assembled';
 if(scale==='face')return {title:'Generic adult face · oral anatomy',subtitle:'Authored facial context · illustrative registration · expert review pending'};
 if(scale==='mouth')return {title:'Oral cavity · selected FDI 36',subtitle:'Authored arches, palate and tongue · no measured tooth registration'};
 if(scale==='tooth')return {title:`FDI 35–37 tooth segment · ${health==='healthy'?'Healthy':'Periodontitis'} · ${form}`,subtitle:'Selected FDI 36 · partial neighbors are authored context · no physical scale'};
 if(scale==='biofilm')return {title:'Biofilm detail · '+PRESETS[health].title,subtitle:'Across attached plaque · schematic observations or separate saved grid'};
 return {title:PRESETS[health].title+' · '+form,subtitle:'Selected distal pocket of FDI 36 · illustrative landmarks · no physical scale'};
}
