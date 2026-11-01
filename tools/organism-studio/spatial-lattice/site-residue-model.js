// Presentation of already qualified correspondences. Evidence qualification stays in the packet reader.
import {endpointLocation,regionResidues} from './interaction-model.js';
export function joinSiteResidues(rows,correspondence){
 const key=r=>JSON.stringify([r.chain,r.number,r.insertion,r.name]),by=new Map();
 for(const row of rows){if(!Number.isInteger(row.index)||row.index<0||!Array.isArray(row.point)||row.point.length!==4||row.point[0]!==row.number||row.point.slice(1).some(v=>!Number.isFinite(v)))throw Error('Source residue coordinates are incomplete; labels held.');const k=key(row);if(by.has(k))throw Error('Ambiguous source residue identity; labels held.');by.set(k,row);}
 const seen=new Set();return correspondence.map(c=>{
  if(!Number.isInteger(c.paper)||!Number.isInteger(c.canonical)||!Number.isInteger(c.author)||typeof c.chain!=='string'||typeof c.insertion!=='string'||typeof c.name!=='string')throw Error('Incomplete numbering correspondence; labels held.');
  const k=JSON.stringify([c.chain,c.author,c.insertion,c.name]),row=by.get(k);
  if(!row||seen.has(k))throw Error('Residue name/chain/author/insertion correspondence differs; labels held.');seen.add(k);
  return {...row,paper:c.paper,canonical:c.canonical,author:c.author,authorLabel:c.chain+':'+c.author+c.insertion,label:c.name[0]+c.name.slice(1).toLowerCase()+c.author+c.insertion+' · '+c.chain};
 });
}
export function siteResidueEntries(data,ledger,location,region){
 if(!region)return [];
 const current=endpointLocation(data,ledger,location.recordId,location.endpointId).regions.find(r=>r.id===region.id);
 if(!current||current.modelProteinId!==region.modelProteinId||JSON.stringify(current.correspondence)!==JSON.stringify(region.correspondence))throw Error('Site selection differs from its qualified source region.');
 return joinSiteResidues(regionResidues(data,ledger,location.recordId,location.endpointId,region.id),current.correspondence).map(r=>({...r,annotationId:current.annotation.id,modelProteinId:current.modelProteinId}));
}
export function siteResiduePage(rows,{query='',page=0,size=24}={}){
 const q=String(query).trim().toLowerCase(),all=rows.filter(r=>!q||[r.name,r.letter,r.fullName,r.chain,r.authorLabel,'paper '+r.paper,'canonical '+r.canonical].some(v=>String(v).toLowerCase().includes(q))),pages=Math.max(1,Math.ceil(all.length/size)),index=Math.min(Math.max(0,page),pages-1);
 return {items:all.slice(index*size,(index+1)*size),total:all.length,page:index,pages,size};
}
export function siteLabelCandidates(rows,selectedIndex=null,maximum=12){
 if(rows.length<=maximum)return [...rows].sort((a,b)=>Number(b.index===selectedIndex)-Number(a.index===selectedIndex));
 const selected=rows.find(r=>r.index===selectedIndex),result=selected?[selected]:[];
 for(let i=0;i<maximum&&result.length<maximum;i++){const row=rows[Math.round(i*(rows.length-1)/(maximum-1))];if(!result.includes(row))result.push(row);}return result;
}
export function residueEvidenceClaim(region){
 const category=region?.category||region?.annotation?.category,status=region?.annotation?.status;
 if(status==='conditional_reference_projection')return {title:'Conditional reference projection',claim:'This residue is projected onto a qualified carrier reference. It is not a measured residue-specific contact or authenticated result for a different strain.'};
 if(category==='assayed_construct')return {title:'Assayed construct membership',claim:'This residue lies within the tested construct. The assay does not identify this individual residue as a contact or binding determinant.'};
 if(category==='measured_reactive_bond')return {title:'Measured reactive-bond residue',claim:'This residue belongs to the experimentally identified reactive bond. This is not a solved protein-pair contact interface.'};
 if(category==='measured_secondary_cleavage')return {title:'Measured secondary-cleavage residue',claim:'This residue belongs to the separately measured secondary cleavage. It is not assigned as an inhibitory site.'};
 if(category==='prediction')return {title:'Predicted region membership',claim:'This is a prediction, not measured interaction or contact-residue evidence.'};
 return {title:'Qualified region membership',claim:'The source supports this region under its recorded qualifications; no residue-specific atomic contact is inferred.'};
}
