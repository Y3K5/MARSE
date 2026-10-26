export type Health = 'healthy'|'periodontitis';
export type Scale = 'pocket'|'tooth'|'mouth'|'face'|'biofilm'|'crevice';
export type Layer = 'structure'|'species'|'environment';
export type Basis = 'evidence'|'model';
export const SOURCES = {
 plaque:{title:'Oral Biofilm Architecture on Natural Teeth · Zijnge et al., 2010',url:'https://pmc.ncbi.nlm.nih.gov/articles/PMC2827546/',method:'FISH observations on seven teeth from four people with advanced periodontitis.',limit:'Selected specimens; static observations, not quantitative abundance or time-resolved succession.'},
 anatomy:{title:'The dimensions of the human dentogingival junction · Vacek et al., 1994',url:'https://pubmed.ncbi.nlm.nih.gov/7928131/',method:'Histomorphometry of 171 surfaces from ten adult cadaver jaws.',limit:'Variable histological dimensions. These cohort observations do not calibrate the geometry shown here.'},
 oxygen:{title:'Oxygen tension in untreated human periodontal pockets',url:'https://pubmed.ncbi.nlm.nih.gov/6592325/',method:'Pocket-base measurements in 111 untreated pockets from 26 participants with advanced disease.',limit:'A base measurement is not a continuous gradient, tissue oxygen measurement or concentration in mol/m³.'},
 fluid:{title:'Flow and albumin content of early pre-inflammatory gingival crevicular fluid',url:'https://doi.org/10.1016/0003-9969(85)90079-2',method:'Fluid sampling during early plaque accumulation and at gingivitis sites.',limit:'No universal flow rate or advanced-pocket calibration is assigned.'}
} as const;
export const TAXA = [
 {id:'ss',name:'Streptococcus sanguinis',short:'S. sanguinis',color:'#7ecfc0',level:'species',region:'unresolved',observation:'This reference does not establish a species-specific subgingival position for S. sanguinis. Streptococcus group observations are not equivalent.'},
 {id:'ac',name:'Actinomyces spp. (group)',short:'Actinomyces spp.',color:'#d6c29a',level:'genus/group',region:'basal',observation:'Actinomyces group cells were observed in the tooth-facing basal portion of attached subgingival biofilm.'},
 {id:'vp',name:'Veillonella parvula',short:'V. parvula',color:'#82a8e4',level:'species',region:'unresolved',observation:'A compatible species-specific location is not established by the selected spatial reference.'},
 {id:'fn',name:'Fusobacterium nucleatum',short:'F. nucleatum',color:'#e7b761',level:'species',region:'intermediate',observation:'F. nucleatum was observed in the intermediate portion of attached biofilm in the selected human specimens.'},
 {id:'pi',name:'Prevotella intermedia',short:'P. intermedia',color:'#bd9ad8',level:'species',region:'microcolonies',observation:'The study describes P. intermedia microcolonies within already formed subgingival biofilms, without a universal pocket-depth location.'},
 {id:'pg',name:'Porphyromonas gingivalis',short:'P. gingivalis',color:'#e67b86',level:'species',region:'microcolonies',observation:'P. gingivalis was described in microcolonies within established subgingival biofilms; this does not establish its position at every site.'},
 {id:'tf',name:'Tannerella forsythia',short:'T. forsythia',color:'#d99764',level:'species',region:'intermediate',observation:'The study reports T. forsythia in the intermediate attached-biofilm region, rather than a universal band at the pocket bottom.'},
 {id:'td',name:'Treponema denticola',short:'T. denticola',color:'#dce779',level:'species',region:'unresolved',observation:'Group-level spirochete observations outside attached biofilm do not establish a T. denticola-specific position. Location remains unresolved.'},
 {id:'pm',name:'Parvimonas micra',short:'P. micra',color:'#9bd0e6',level:'species',region:'microcolonies',observation:'P. micra is discussed among taxa forming microcolonies in established biofilm; exact site-specific positions and abundance remain unresolved.'}
] as const;
export const PRESETS = {
 healthy:{margin:.35,attachment:-1.65,crest:-3.5,title:'Healthy sulcus'},
 periodontitis:{margin:-.15,attachment:-6.2,crest:-8.0,title:'Periodontitis pocket'}
} as const;
export function observation(id:string,health:Health){
 const taxon=TAXA.find(t=>t.id===id);
 if(!taxon) return {resolved:false,region:'unresolved',text:'Select a taxon to inspect its source and identification level.'};
 if(health==='healthy') return {resolved:false,region:'unresolved',text:'Location unresolved for healthy anatomy. The selected spatial study concerns advanced periodontitis and is not transferred to a healthy sulcus.'};
 return {resolved:taxon.region!=='unresolved',region:taxon.region,text:taxon.observation};
}
export interface SavedCase {id:Health;frame_index:number;model_time:number;nx:number;ny:number;biomass:number[][];}
export interface SavedInput {format:string;host:'human';taxa:Array<{id:string;name:string;short:string;color:string}>;cases:SavedCase[];source:{file:string;sha256:string;precision:string};units:string;coordinates:{order:string;x:string;y:string;physical_registration:null};}
export function validateSaved(data:SavedInput):SavedInput {
 if(data.format!=='marse.pocket-explorer-fields/1'||data.host!=='human'||!/^[a-f0-9]{64}$/.test(data.source?.sha256))throw Error('Incompatible source identity');
 if(data.units!=='arbitrary model density'||data.coordinates?.order!=='row-major y*nx+x'||data.coordinates.x!=='Across the assumed model channel'||data.coordinates.y!=='Axial model depth; y=0 at the entrance, increasing apically'||data.coordinates.physical_registration!==null||!data.source.precision)throw Error('Incompatible units, axes or precision');
 if(data.taxa.length!==TAXA.length||data.taxa.some((t,i)=>t.id!==TAXA[i].id||t.color!==TAXA[i].color))throw Error('Incompatible taxon identities');
 if(data.cases.length!==2||new Set(data.cases.map(c=>c.id)).size!==2)throw Error('Missing source cases');
 for(const c of data.cases){if(!['healthy','periodontitis'].includes(c.id)||c.nx!==32||c.ny!==72||c.biomass.length!==9||!Number.isInteger(c.frame_index)||!Number.isFinite(c.model_time)||c.biomass.some(a=>a.length!==c.nx*c.ny||a.some(v=>!Number.isFinite(v)||v<0)))throw Error('Incompatible density grid');}
 return data;
}
export function gridRGBA(run:SavedCase,id:string){
 const index=TAXA.findIndex(t=>t.id===id);if(index<0&&id!=='all')throw Error('Unknown taxon');
 const rgb=TAXA.map(t=>t.color.match(/[a-f0-9]{2}/gi)!.map(v=>parseInt(v,16)));
 const values=Array.from({length:run.nx*run.ny},(_,i)=>index<0?run.biomass.reduce((s,a)=>s+a[i],0):run.biomass[index][i]);
 const max=Math.max(...values)||1,bytes=new Uint8ClampedArray(run.nx*run.ny*4);
 for(let i=0;i<values.length;i++){
  let k=index;if(k<0){k=0;for(let j=1;j<9;j++)if(run.biomass[j][i]>run.biomass[k][i])k=j;}
  bytes.set([...rgb[k],values[i]>0?Math.round(85+170*Math.sqrt(values[i]/max)):0],i*4);
 }
 return bytes;
}
