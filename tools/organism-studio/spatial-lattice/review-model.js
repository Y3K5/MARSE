// Navigation over the validated retained packet. Never creates an interaction edge.
const caches=new WeakMap();
const namedPairLevels=new Set(['direct assay','direct biochemical']);
// Evidence classification only; this does not assign a retained reference, strain or contact site.
export function isSupportedProteinPair(data,record){
 return namedPairLevels.has(record.level)&&[record.a,record.b].every(id=>data.proteins.some(p=>p.id===id));
}
const tokens=s=>String(s||'').toLowerCase().split(/[^a-z0-9]+/).filter(Boolean);
function classify(data,r,other){
 if(data.contexts.some(p=>p.id===other.id))return 'unresolved_partner';
 if(r.level==='candidate / unresolved identity')return 'candidate';
 return isSupportedProteinPair(data,r)?'supported_pair':'source_context';
}
// Consume producer-validated full-reference correspondences, never labels or homologs.
export function sourceReferences(data,id){
 const p=data.proteins.find(v=>v.id===id);if(!p?.subtractive_record)return [];
 return (data.subtractive?.reference_links||[]).filter(l=>l.a===id&&l.kind==='exact_reference_identity'&&typeof l.sequence_sha256==='string'&&!!l.sequence_sha256&&l.sequence_sha256===p.subtractive_record.sequence_sha256&&data.proteins.some(v=>v.id===l.b&&!v.subtractive_record&&v.species===p.species));
}
export function reviewAudit(data){
 if(caches.has(data))return caches.get(data);
 const keys=new Set(data.subtractive?.shortlist_keys||[]),by=new Map([...data.proteins,...data.contexts].map(p=>[p.id,p]));
 const entries=data.proteins.filter(p=>keys.has(p.subtractive_record?.key)).map(p=>{
  const refs=sourceReferences(data,p.id);
  const mapped=[];for(const link of refs)for(const r of data.relationships){const otherId=r.a===link.b?r.b:r.b===link.b?r.a:null,other=by.get(otherId);if(other&&other.species!==p.species)mapped.push({recordId:r.id,sourceEndpoint:link.b,otherEndpoint:other.id,referenceLink:link.id,category:classify(data,r,other)});}
  // Same labels offer source navigation only. No homolog/strain correspondence is asserted.
  const labelTokens=tokens(p.name),sourceIds=new Set(data.subtractive.source_context_ids),holds=[];
  for(const source of data.proteins){if(!sourceIds.has(source.id)||source.species!==p.species||refs.some(l=>l.b===source.id))continue;const nameTokens=tokens(source.name);if(nameTokens.length!==1||nameTokens[0].length<4||!labelTokens.includes(nameTokens[0]))continue;
   const records=data.relationships.filter(r=>[r.a,r.b].includes(source.id)&&by.get(r.a)?.species!==by.get(r.b)?.species);if(!records.length)continue;
   holds.push({sourceEndpoint:source.id,recordIds:records.map(r=>r.id),category:'label_only_mapping_hold',reason:'Similar label only. The retained reference does not pass the exact full-sequence endpoint gate. Source strain/material and interaction evidence must not transfer to this review protein.',sourceAccession:source.accession,sourceStrain:source.reference?.organism?.scientificName||'Unknown source strain',sourceLength:source.reference?.length,queryLength:p.subtractive_record.length});
  }
  return {proteinId:p.id,key:p.subtractive_record.key,name:p.name,mapped,holds,status:mapped.some(r=>r.category==='supported_pair')?'supported_pair':mapped.length?'source_context':holds.length?'mapping_hold':'unresolved'};
 });
 const summary={total:entries.length,retainedFolds:entries.filter(e=>data.shape_traces?.[e.proteinId]).length,coordinateGaps:entries.filter(e=>!data.shape_traces?.[e.proteinId]).length,supportedProteins:entries.filter(e=>e.status==='supported_pair').length,supportedRecords:[...new Set(entries.flatMap(e=>e.mapped.filter(r=>r.category==='supported_pair').map(r=>r.recordId)))].length,mappedContextProteins:entries.filter(e=>e.status==='source_context').length,heldProteins:entries.filter(e=>e.status==='mapping_hold').length,unresolvedProteins:entries.filter(e=>e.status==='unresolved').length,entries};caches.set(data,summary);return summary;
}
export function reviewEntry(data,id){return reviewAudit(data).entries.find(e=>e.proteinId===id)||null;}
export function reviewStyle(entry){return entry?.status==='supported_pair'?{key:'supported',color:'#78dcc5',label:'Review · source-supported pair'}:entry?.status==='mapping_hold'?{key:'hold',color:'#f3ae78',label:'Review · source mapping hold'}:{key:'review',color:'#f1d47b',label:'Review · connection unresolved'};}
