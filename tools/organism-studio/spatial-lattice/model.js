import {displayLimits} from './exploration-model.js';
import {validatePairFocus} from './pair-model.js';
import {validateNeighborhood,validateNeighborhoodSelection} from './neighborhood-model.js';
import {validateNavigationProtein} from './navigation-model.js';
export const palette={},basePositions={unknown:[-35,35,30]};
export function configureSpecies(data){const colors=['#e6a08c','#e0c07d','#88cbbd','#afc3ef','#d5a3d9'];for(const [i,s] of data.species.entries()){palette[s.id]=/^#[0-9a-f]{6}$/i.test(s.color||'')?s.color:colors[i%colors.length];const angle=i/Math.max(1,data.species.length)*Math.PI*2;basePositions[s.id]=[Math.cos(angle)*40,(i%3-1)*12,Math.sin(angle)*40];}palette.unknown='#b5bdc4';}
export function positionFor(id,spacing=1){return basePositions[id]?.map(x=>x*spacing)||null;}
export function sphereSlot(i,n,r=14){const y=1-2*(i+.5)/n,a=i*Math.PI*(3-Math.sqrt(5));return [Math.cos(a)*Math.sqrt(1-y*y)*r,y*r,Math.sin(a)*Math.sqrt(1-y*y)*r];}
export function traceFrame(shape){
 const rows=shape.chains.flatMap(c=>c.ca);if(!rows.length)throw Error('No retained source coordinates.');
 const center=[0,0,0];for(const r of rows)for(let j=0;j<3;j++)center[j]+=r[j+1]/rows.length;
 let radius=0;for(const r of rows)radius=Math.max(radius,Math.hypot(...center.map((v,j)=>r[j+1]-v)));
 const segments=[];for(const c of shape.chains){let part=[];for(const r of c.ca){const last=part.at(-1);if(last&&(r[0]-last[0]<=0||r[0]-last[0]>2||Math.hypot(...[1,2,3].map(j=>r[j]-last[j]))>6)){if(part.length>1)segments.push({chain:c.chain,rows:part});part=[];}part.push(r);}if(part.length>1)segments.push({chain:c.chain,rows:part});}
 return {center,radius,segments,count:rows.length};
}
export function niceStep(distance){const value=Math.max(.01,distance/8),power=10**Math.floor(Math.log10(value)),fraction=value/power;return (fraction<=1?1:fraction<=2?2:fraction<=5?5:10)*power;}
export function assertMap(d){
 const safe=v=>typeof v==='string'&&/^[a-zA-Z0-9_.:|/-]+$/.test(v)&&!v.includes('..')&&!['__proto__','constructor','prototype'].includes(v);
 if(d.schema!=='grouped-organism-workspace-v1'||d.native_placements_admitted!==0||d.candidate_interaction_transfer!==false)throw Error('Unsupported source-map contract.');
 for(const key of ['species','proteins','contexts','systems','relationships','sources','placements'])if(!Array.isArray(d[key]))throw Error('Missing source-map collection: '+key);
 for(const key of ['pairs','shape_traces','workspace_bridge','study_contexts'])if(!d[key]||typeof d[key]!=='object'||Array.isArray(d[key]))throw Error('Missing source-map object: '+key);
 const all=[...d.proteins,...d.contexts],ids=new Set(all.map(x=>x.id));if(ids.size!==all.length||all.some(p=>!safe(p.id)||!safe(p.species)))throw Error('Invalid or duplicate endpoint identity.');
 for(const p of d.proteins){if(p.reference?.sequence!=null&&typeof p.reference.sequence!=='string'||p.reference?.length!=null&&(!Number.isInteger(p.reference.length)||p.reference.length<1)||p.subtractive_record&&(!Number.isInteger(p.subtractive_record.length)||p.subtractive_record.length<1))throw Error('Invalid reference sequence/length.');}
 const taxa=new Set(d.species.map(s=>s.id));if(taxa.size!==d.species.length||d.species.some(s=>!safe(s.id))||all.some(p=>!taxa.has(p.species)&&p.species!=='unknown'))throw Error('Unresolved species reference.');
 const records=new Set(d.relationships.map(r=>r.id));if(records.size!==d.relationships.length)throw Error('Duplicate relationship identity.');
 for(const r of d.relationships)if(!safe(r.id)||!ids.has(r.a)||!ids.has(r.b))throw Error('Relationship endpoint is missing.');
 for(const g of d.systems)if(!safe(g.id)||!Array.isArray(g.members)||g.members.some(id=>!ids.has(id)))throw Error('System membership differs.');
 for(const pair of Object.values(d.pairs))if(!safe(pair.id)||pair.record_ids.some(id=>!records.has(id)))throw Error('Pair record is missing.');
 for(const source of d.sources){const u=new URL(source.url);if(!['http:','https:'].includes(u.protocol)||u.username||u.password)throw Error('Unsafe source URL.');}
 for(const [id,shape]of Object.entries(d.shape_traces)){if(!ids.has(id)||!Array.isArray(shape.chains)||!shape.chains.length||!Number.isFinite(traceFrame(shape).radius)||traceFrame(shape).radius<=0)throw Error('Invalid retained shape.');}
 if(d.subtractive){const sub=d.subtractive,rows=d.proteins.filter(p=>p.subtractive_record),keys=rows.map(p=>p.subtractive_record.key);if(sub.schema!=='retained-subtractive-marse-v1'||rows.length!==sub.record_count||new Set(keys).size!==keys.length||!Array.isArray(sub.shortlist_keys)||sub.shortlist_keys.some(k=>!keys.includes(k))||!Array.isArray(sub.reference_links)||sub.reference_links.some(l=>l.kind!=='exact_reference_identity'||!ids.has(l.a)||!ids.has(l.b)))throw Error('Result or reference-link membership differs.');}
 return d;
}
export function validateView(s,data){
 const finite=v=>Array.isArray(v)&&v.length===3&&v.every(x=>Number.isFinite(x)&&Math.abs(x)<1e6);
 if(s?.schema!=='marse-lattice-view'||s.version!==1||!['microbes','desk','anatomy'].includes(s.preset)||!['overview','network','organism','protein','pen','cell'].includes(s.context?.kind))throw Error('Unsupported lattice view.');
 if(s.sourceSHA256!==data.sha256)throw Error('Source snapshot does not match this workspace.');
 const kinds={microbes:['overview','network','organism','protein'],desk:['overview','pen'],anatomy:['overview','cell']};if(!kinds[s.preset].includes(s.context.kind))throw Error('Reference frame does not belong to this scene.');
 if(s.context.kind==='organism'&&!data.map.species.some(p=>p.id===s.context.id))throw Error('Unknown organism frame.');
 if(s.context.kind==='protein'&&!Object.hasOwn(data.map.shape_traces,s.context.id))throw Error('No retained coordinates for this protein.');
 if(s.context.returnTo!==undefined&&s.context.returnTo!=='network')throw Error('Invalid source-frame return.');
 if(s.networkDisplay&&(!['all','across','direct'].includes(s.networkDisplay.scope)||['candidates','bodies','unlinked'].some(k=>typeof s.networkDisplay[k]!=='boolean')))throw Error('Invalid network options.');
 if(s.networkDisplay?.subtractive){const v=s.networkDisplay.subtractive;if(!data.map.subtractive||!['forms','shortlist','catalogue'].includes(v.cohort)||!['all',...data.map.subtractive.species.map(s=>s.id)].includes(v.organism)||!['all',...data.map.subtractive.groups.map(g=>g.id)].includes(v.group)||typeof v.query!=='string'||v.query.length>2000||!Number.isInteger(v.page)||v.page<0||v.page>1000||typeof v.sourceContext!=='boolean'||typeof v.referenceLinks!=='boolean'||v.reviewHighlights!==undefined&&typeof v.reviewHighlights!=='boolean')throw Error('Invalid saved subtractive filters.');}
 if(s.networkDisplay?.pairFocus){validatePairFocus(data.map,s.networkDisplay.pairFocus);if(s.context.kind!=='network'||s.preset!=='microbes')throw Error('Large-pair display requires the network frame.');const v=s.networkDisplay.subtractive;if(data.map.subtractive&&(!v||v.cohort!=='forms'||v.organism!=='all'||!v.sourceContext))throw Error('Large-pair display requires the complete source context.');}
 if(s.networkDisplay?.neighborhood!==undefined){validateNeighborhood(data.map,s.networkDisplay.neighborhood);if(s.networkDisplay.pairFocus||s.context.kind!=='network'||s.preset!=='microbes')throw Error('Protein neighborhood requires its own network frame.');}
 if(s.networkDisplay?.navigationProtein!==undefined){if(s.networkDisplay.pairFocus||!['network','protein'].includes(s.context.kind)||s.context.kind==='protein'&&s.context.id!==s.networkDisplay.navigationProtein)throw Error('Saved selection conflicts with the visible frame.');if(s.networkDisplay.neighborhood)validateNeighborhoodSelection(data.map,s.networkDisplay.neighborhood,s.networkDisplay.navigationProtein);else validateNavigationProtein(data.map,s.networkDisplay.subtractive||{},s.networkDisplay.navigationProtein);}
 if(s.networkDisplay?.overviewCamera){const c=s.networkDisplay.overviewCamera;if(!finite(c.position)||!finite(c.target)||Math.hypot(...c.position.map((v,i)=>v-c.target[i]))<.2)throw Error('Invalid overview camera.');}
 if(s.networkDisplay?.layout!==undefined&&!['cluster','layers','orbit'].includes(s.networkDisplay.layout))throw Error('Unknown microbial neighborhood layout.');
 if(s.networkDisplay?.foldSize!==undefined&&(!Number.isFinite(s.networkDisplay.foldSize)||s.networkDisplay.foldSize<displayLimits.foldMin||s.networkDisplay.foldSize>displayLimits.foldMax))throw Error('Invalid fold display size.');
 if(s.proteinDisplay&&(!['path','amino','chain'].includes(s.proteinDisplay.colors)||typeof s.proteinDisplay.forms!=='boolean'||typeof s.proteinDisplay.beads!=='boolean'))throw Error('Invalid protein display options.');
 if(!finite(s.camera?.position)||!finite(s.camera?.target)||Math.hypot(...s.camera.position.map((v,i)=>v-s.camera.target[i]))<.2)throw Error('Invalid camera.');
 if(!Number.isFinite(s.spacing)||s.spacing<displayLimits.spacingMin||s.spacing>displayLimits.spacingMax||!Number.isFinite(s.explode)||s.explode<0||s.explode>1||!['high','balanced'].includes(s.quality)||typeof s.grid!=='boolean'||typeof s.labels!=='boolean'||typeof s.cutaway!=='boolean')throw Error('Invalid view options.');return s;
}
