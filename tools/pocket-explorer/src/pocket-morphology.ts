import {TAXA} from './pocket-science';
export type CellForm='cocci'|'branching'|'fusiform'|'short-rod'|'long-rod'|'spirochete';
export const MORPHOLOGY:Record<string,{form:CellForm;label:string;observation:string;source:string;context:string}>={
 ss:{form:'cocci',label:'Cocci · illustrative chain',observation:'Coccal chains and clusters; chain length varies. The number of cells drawn is arbitrary.',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC3919539/',context:'Phase-contrast study of S. sanguinis SK36; morphology is strain- and condition-dependent.'},
 ac:{form:'branching',label:'Rods / branching filaments',observation:'An illustrative branching form for the Actinomyces group. Members can also appear as short rods; this is not a species identification.',source:'https://www.sgmjournals.org/mic/content/155/7/2116',context:'Human initial-plaque microscopy of A. naeslundii provides an example, not a shape for every group member.'},
 vp:{form:'cocci',label:'Small cocci',observation:'Rounded coccal form. The pair shown is a display choice, not a measured association.',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC3035260/',context:'V. parvula type-strain genome report with morphology description.'},
 fn:{form:'fusiform',label:'Elongated fusiform rod',observation:'An elongated, tapered cell body. Shape alone cannot identify this taxon or establish its partners.',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC2818713/',context:'TEM of a laboratory oral consortium including F. nucleatum subsp. vincentii; no laboratory organization is transferred into the saved grid.'},
 pi:{form:'short-rod',label:'Rod-shaped cell',observation:'Rounded rod body. No membrane-protein sites or extracellular-vesicle counts are drawn.',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC10897501/',context:'P. intermedia cryo-electron tomography; only qualitative body shape is represented.'},
 pg:{form:'short-rod',label:'Short rod / ovoid rod',observation:'A short, plump rod rather than a perfect sphere. Surface relief is artistic, not a measured envelope.',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC10440386/',context:'Untreated-control SEM observations in a laboratory study; no treatment effects or experimental conditions are modeled.'},
 tf:{form:'long-rod',label:'Elongated rod',observation:'An elongated cell body. Natural-tooth microscopy also documents filamentous Tannerella forms. No crystalline S-layer spacing is claimed.',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC3354324/',context:'T. forsythia ultrastructure study. This surface illustration is not a molecular reconstruction.'},
 td:{form:'spirochete',label:'Irregular twisted spirochete',observation:'A thin, irregularly twisted cell, with planar and helical regions possible. A perfectly regular spring would conceal that variability.',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC178876/',context:'Microscopy of T. denticola cell morphology; it does not establish a species-specific pocket location.'},
 pm:{form:'cocci',label:'Small cocci',observation:'Coccal bodies; the number and arrangement drawn are illustrative.',source:'https://www.frontiersin.org/journals/public-health/articles/10.3389/fpubh.2022.994279/full',context:'Clinical isolate microscopy from an extraoral case; used for morphology only, not oral colonization evidence.'}
};

// Authored vector surface cues. Colors identify taxa; they are not microscopy
// colors. Shapes are normalized separately to fit, not drawn at relative sizes.
export function morphologySVG(id:string,compact=false){
 const t=TAXA.find(t=>t.id===id);if(!t)throw Error('Unknown morphology identity');
 const m=MORPHOLOGY[id],key='cell-'+id;
 const circle=(x:number,y:number,r=15)=>`<circle cx="${x}" cy="${y}" r="${r}" fill="url(#${key})" stroke="${t.color}" stroke-opacity=".38" stroke-width=".8"/>`;
 let shape='';
 if(m.form==='cocci')shape=id==='ss'?[-58,-29,0,29,58].map((x,i)=>circle(150+x,70+Math.sin(i*.8)*6)).join(''):id==='vp'?circle(133,74,22)+circle(175,64,22):circle(135,72,24)+circle(184,80,18);
 if(m.form==='short-rod')shape='<rect x="84" y="48" width="132" height="49" rx="24.5" fill="url(#'+key+')" stroke="'+t.color+'" stroke-width=".7"/>';
 if(m.form==='long-rod')shape='<path d="M37 74Q64 53 144 55Q227 53 263 73Q226 92 146 91Q64 93 37 74Z" fill="url(#'+key+')"/>';
 if(m.form==='fusiform')shape='<path d="M23 73Q72 72 94 61Q130 49 173 57Q230 71 277 74Q230 76 191 86Q149 99 104 87Q66 77 23 73Z" fill="url(#'+key+')"/>';
 if(m.form==='branching')shape='<path d="M55 108Q89 73 143 73Q190 71 247 40M116 75Q104 48 85 30M181 69Q194 91 222 106" stroke="url(#'+key+')" stroke-width="13" stroke-linecap="round" fill="none"/>';
 if(m.form==='spirochete')shape='<path d="M27 75C39 30 50 117 66 74S91 25 102 66S118 108 138 73S160 28 173 79S194 121 209 66S240 35 253 77S269 93 277 67" fill="none" stroke="url(#'+key+')" stroke-width="5.5" stroke-linecap="round"/>';
 return `<svg class="cell-form${compact?' compact':''}" viewBox="0 0 300 146" role="img" aria-label="${t.name}: ${m.label}; schematic morphology, normalized display size"><defs><linearGradient id="${key}" x1="0" y1="0" x2="0.2" y2="1"><stop stop-color="#ecf1df"/><stop offset=".25" stop-color="${t.color}"/><stop offset=".62" stop-color="${t.color}"/><stop offset="1" stop-color="#203137"/></linearGradient></defs>${shape}</svg>`;
}
export function morphologyPanel(id:string){
 if(id==='all')return `<section class="morphology-panel"><h3>Species form library</h3><p class="morphology-limit">Qualitative forms · normalized to fit · no relative size or cell counts</p><div class="cell-gallery">${TAXA.map(t=>`<button data-morphology="${t.id}" aria-label="View ${t.name} structure">${morphologySVG(t.id,true)}<span>${t.short}</span><small>${MORPHOLOGY[t.id].label}</small></button>`).join('')}</div><p class="morphology-limit">Select a form to see its microscopy source. These cells are not placed in the map.</p></section>`;
 const t=TAXA.find(t=>t.id===id);if(!t)throw Error('Unknown morphology identity');const m=MORPHOLOGY[id];
 return `<section class="morphology-panel"><div class="eyebrow">Species structure · qualitative reference</div><h3>${t.name}</h3>${morphologySVG(id)}<p class="cell-form-label">${m.label}</p><p>${m.observation}</p><p class="morphology-limit">${m.context}</p><a href="${m.source}" target="_blank" rel="noopener noreferrer">Morphology evidence ↗</a><p class="morphology-limit">Artistic color and relief · no physical scale, counted cells, molecular binding sites or inferred contacts. Saved density does not determine cell shape or location.</p></section>`;
}
