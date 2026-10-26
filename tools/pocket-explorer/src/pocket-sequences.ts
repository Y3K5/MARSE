import {REFERENCES,type ReferenceKey} from './pocket-reference';
export type Sequence='assembly'|'response';
export interface Stage {title:string;status:'Observed'|'Inferred'|'Illustrative';text:string;source:ReferenceKey;}
export const SEQUENCES:Record<Sequence,Stage[]>={
 assembly:[
  {title:'Acquired pellicle',status:'Illustrative',text:'An acquired salivary film conditions the tooth surface. Its thickness, formation rate and molecular composition are not modeled here.',source:'supra'},
  {title:'Surface attachment',status:'Inferred',text:'Anonymous bacterial shapes attach to the conditioned surface. This explanatory transition does not assign species, binding strengths or measured contact sites.',source:'supra'},
  {title:'Local biomass accumulation',status:'Inferred',text:'Local growth and recruitment can contribute to an aggregate. The visible glyphs illustrate a process; their count, expansion rate and relative sizes are arbitrary.',source:'supra'},
  {title:'Coaggregation and matrix',status:'Illustrative',text:'Schematic neighboring cells and an extracellular matrix illustrate assembly. No specific adhesion pair, metabolic flux or matrix yield is calculated.',source:'supra'},
  {title:'Organized community',status:'Observed',text:'Spatial microscopy documents organized mature communities. The anonymous animation is illustrative; the separate Evidence view retains specimen- and taxon-specific observations.',source:'plaque'},
  {title:'Turnover and dispersal',status:'Illustrative',text:'Community turnover is illustrated by a few glyphs moving away. Detachment frequency, fluid flow and microbial survival are unresolved.',source:'supra'}
 ],
 response:[
  {title:'Tooth–gingiva barrier',status:'Illustrative',text:'The epithelial boundary and connective tissues are separate compartments. This storyboard does not move the anatomical attachment or change the chosen preset.',source:'attachment'},
  {title:'Plaque-associated gingival response',status:'Observed',text:'Experimental human gingivitis provides evidence of a gingival response to plaque accumulation. Its findings do not provide a clock for periodontitis or bone loss.',source:'gingivitis'},
  {title:'Inflammatory signaling',status:'Illustrative',text:'A tissue-side glow represents a host inflammatory response. It is not a cytokine concentration, immune-cell count or receptor-binding result.',source:'gingivitis'},
  {title:'Altered local ecology',status:'Inferred',text:'Host responses and microbial ecology can interact. Oxygen, pH, nutrient supply and perfusion at this illustrated site remain unknown.',source:'plaque'},
  {title:'Attachment loss · possible disease outcome',status:'Illustrative',text:'Attachment loss is shown as a conceptual landmark in a separate storyboard. It is not an inevitable outcome of these earlier stages and is not a prediction for this tooth.',source:'attachment'},
  {title:'Bone loss · separate endpoint',status:'Illustrative',text:'A separate conceptual bone outline explains loss of supporting tissue. No osteoclast activity, remodeling rate or clinical progression is calculated.',source:'boneDefects'}
 ]
};
export function stageAt(sequence:Sequence,progress:number){
 if(!Number.isFinite(progress))throw Error('Invalid storyboard position');
 const index=Math.min(5,Math.max(0,Math.floor(progress)));return {...SEQUENCES[sequence][index],index};
}
// Deterministic visual choreography, explicitly not a biological simulation.
export function sequenceSVG(sequence:Sequence,progress:number){
 const t=Math.max(0,Math.min(5,progress)),k=Math.floor(t),stage=stageAt(sequence,t);
 const fill=sequence==='response'?'#9b5c50':'#727c59';
 let glyphs='';
 for(let i=0;i<30;i++){
  const threshold=1+(i%7)*.35,alpha=Math.min(1,Math.max(0,(t-threshold)*2));
  const x=126+(i%6)*53+Math.sin(i*2.3)*12,y=276-Math.floor(i/6)*30;
  const drift=t>4.5&&i%9===0?(t-4.5)*105:0;
  glyphs+=`<g opacity="${alpha*(drift?1-drift/90:1)}" transform="translate(${x+drift},${y-drift}) rotate(${i*43%135})"><rect x="-12" y="-5" width="24" height="10" rx="5" fill="${i%3===0?'#c5d5c1':i%3===1?'#c8b384':'#d4a89b'}"/><path d="M-7 -2H7" stroke="#fff" stroke-opacity=".18"/></g>`;
 }
 const matrix=Math.max(0,Math.min(.65,(t-2.2)*.25));
 const response=sequence==='response'?`<path d="M495 90Q466 195 490 290L598 294V87Z" fill="#9a5a54"/><path d="M491 104Q465 189 490 278" stroke="#edb2a2" stroke-width="6" fill="none"/><ellipse cx="549" cy="190" rx="52" ry="99" fill="#ef8973" opacity="${Math.min(.4,t*.1)}"/><path d="M495 ${t>=4?326:295}H592V344H495" fill="#dfcbae"/><path d="M493 ${t>=4?325:294}H530" stroke="#bfe2cf" stroke-width="2" stroke-dasharray="4 4"/>`:`<path d="M493 91Q466 191 489 291" stroke="#c1938b" stroke-width="9" fill="none"/>`;
 return `<svg class="sequence-art" viewBox="0 0 640 370" role="img" aria-label="${stage.title}; anonymous illustrative storyboard, no measured cell positions or elapsed time"><defs><linearGradient id="seq-tooth" x2="1" y2=".2"><stop stop-color="#9f927b"/><stop offset=".7" stop-color="#f0e5c9"/><stop offset="1" stop-color="#d8c3a0"/></linearGradient><radialGradient id="seq-matrix"><stop stop-color="${fill}"/><stop offset="1" stop-color="#243a3c"/></radialGradient></defs><rect width="640" height="370" fill="#101e23"/><path d="M60 70Q90 165 72 308H105Q128 205 112 68Z" fill="url(#seq-tooth)"/><path d="M110 91Q126 192 114 294" stroke="#c8b98b" stroke-width="${3+t*.6}" fill="none"/><path d="M115 91Q180 76 272 108T445 101Q465 198 445 285Q293 319 119 295Z" fill="url(#seq-matrix)" opacity="${matrix}"/>${glyphs}${response}<text x="58" y="40" fill="#dbc9a6" font-size="15">Tooth surface</text><text x="504" y="40" fill="#c89e99" font-size="15">${sequence==='response'?'Host tissue':'Fluid-side boundary'}</text><text x="316" y="345" text-anchor="middle" fill="#a6bcb3" font-size="12">Anonymous shapes · no cell counts, physical scale or elapsed-time calibration</text></svg>`;
}
export function sourceForStage(sequence:Sequence,progress:number){return REFERENCES[stageAt(sequence,progress).source];}
