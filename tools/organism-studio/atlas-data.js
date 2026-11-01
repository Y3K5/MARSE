// Viewing dimensions are explicitly illustrative. Molecular files retain Å coordinates.
export const referenceTaxa=[
 {id:'ecoli_k12',name:'Escherichia coli K-12',short:'E. coli K-12',color:'#74b8c5',shape:'rod',set:'Catalogue reference'},
 {id:'bsubtilis',name:'Bacillus subtilis',short:'B. subtilis',color:'#c9b591',shape:'rod',set:'Catalogue reference'},
 {id:'scerevisiae',name:'Saccharomyces cerevisiae',short:'S. cerevisiae',color:'#bb9acd',shape:'yeast',set:'Catalogue reference'}];
export const profiles={
 td:{envelope:'diderm',length:8,width:.28,summary:'Spirochete · internal periplasmic flagella · separate cytoplasmic filament ribbon',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC2519799/'},
 pg:{envelope:'diderm',length:1.2,width:.6,summary:'Short rod · external fimbriae · variable capsule',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC10929438/'},
 tf:{envelope:'diderm',length:3,width:.45,summary:'Fusiform rod · glycoprotein S-layer outside outer membrane',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC3354324/'},
 ss:{envelope:'monoderm',length:.8,width:.8,summary:'Coccal chains · cell-wall-anchored adhesins',source:'https://www.rcsb.org/structure/5EQ2'},
 ac:{envelope:'monoderm',length:3,width:.5,summary:'Filament-form Actinomyces group · A. oris supplies the molecular example',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC4571867/'},
 vp:{envelope:'diderm',length:.5,width:.5,summary:'Diderm coccus · OmpM supports outer-membrane attachment and transport',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC5491611/'},
 fn:{envelope:'diderm',length:5,width:.5,summary:'Fusiform rod · FadA adhesion system',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC2635053/'},
 pi:{envelope:'diderm',length:1.5,width:.6,summary:'Anaerobic rod · secreted interpain A protease example',source:'https://www.rcsb.org/structure/3BBA'},
 pm:{envelope:'monoderm',length:.6,width:.6,summary:'Coccus · surface-associated TmpC evidence; coordinate acquisition unresolved',source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC13080454/'},
 ecoli_k12:{envelope:'diderm',length:2,width:.8,summary:'Catalogue K-12 reference · OmpF porin and FhuA receptor',source:'https://journals.asm.org/doi/10.1128/jb.185.20.6112-6118.2003'},
 bsubtilis:{envelope:'monoderm',length:2.4,width:.9,summary:'Catalogue reference · thick cell wall · extracellular TasA scaffold',source:'https://www.nature.com/articles/s41467-022-34700-z'},
 scerevisiae:{envelope:'yeast',length:5,width:5,summary:'Catalogue reference yeast · plasma membrane, cell wall and internal microtubules',source:'https://www.rcsb.org/structure/5W3H'}};
export const extraProteins={
 ss:[['SrpA binding region','Cell-wall-associated adhesin','Experimental binding-region construct; full serine-rich stalk and cell-wall linkage absent.','EXT_5EQ2','surface','5EQ2']],
 ac:[['FimP pilin domains','Pilus subunit','A. oris T14V mutated soluble-domain construct. Group-level Actinomyces display is not a strain match.','EXT_3UXF','surface','3UXF']],
 vp:[['OmpM1','Outer-membrane porin / tether','Trimeric nutrient channel with peptidoglycan-associated tethering context; not all domains are resolved.','EXT_8BYS','surface','8BYS']],
 fn:[['FadA adhesin','Surface-associated adhesion system','Mature soluble monomer has an antiparallel helical fold. Deposited monomer is not a full native filament.','EXT_3ETW','surface','3ETW']],
 pi:[['Interpain A','Extracellular protease','Active protease-domain structure; precursor and secretion context absent.','EXT_3BBA','surface','3BBA']],
 pm:[['TmpC','Reported surface-associated protein','UniProt A0A0B4S1T1 matches P. micra. AlphaFold lists a model but its download was refused; no substitute species model is shown.','PM_TMP','uncertain']],
 ecoli_k12:[['OmpF porin','Outer-membrane channel','OPM-supplied trimer; individual chains are β-barrels. Not a ligand-specific receptor.','EXT_2OMF','surface','2OMF'],['FhuA receptor','Outer-membrane receptor','Ferrichrome receptor with an OPM membrane frame. Retained LPS is not a reconstructed leaflet.','EXT_2FCP','surface','2FCP']],
 bsubtilis:[['KtrAB transporter','Cytoplasmic-membrane complex','KtrB membrane components and intracellular KtrA regulator; OPM-selected complex, experimental mutations.','EXT_4J7C','internal','4J7C'],['OpuAC binding domain','Peripheral lipoprotein domain','Isolated extracellular binding domain; native anchor/transporter absent.','EXT_2B4L','surface','2B4L'],['TasA monomer','Extracellular scaffold protein','Soluble chain A; the monomer is not a filament.','EXT_5OF1','surface','5OF1'],['TasA biofilm filament','Extracellular polymer','Experimental assembled filament uses donor-strand exchange.','EXT_8AUR','surface','8AUR']],
 scerevisiae:[['Ste2 receptor complex','Plasma membrane beneath cell wall','Engineered receptor–G-protein–pheromone complex; intracellular partners are not exposed.','EXT_7AD3','internal','7AD3'],['Pma1 proton pump','Plasma membrane beneath cell wall','Autoinhibited hexamer; cytosolic domains are retained.','EXT_7VH5','internal','7VH5'],['α/β-tubulin','Internal microtubule subunit','Experimental α/β dimer from stabilized yeast microtubules. Distinct from bacterial flagellin and CfpA.','EXT_5W3H','internal','5W3H']]
};
export const tdMachinery=[
 ['Flagellar filament segment','Periplasmic polymer','42-chain experimental segment with several sheath-protein types and flagellin. The real assembly is more complex than a uniform FlaA coating.','EXT_9PIN','internal','9PIN'],
 ['CfpA ribbon','Cytoplasmic filaments','Filament ribbon runs inside the inner membrane, distinct from periplasmic flagella. Atomic model not assigned here.','TD_CFP','internal'],
 ['Flagellar motor + hook','Inner-membrane motor','Tomography supports polar basal bodies and periplasmic hooks. Motor component shapes are conceptual, not a fitted atomic motor.','TD_MOTOR','internal'],
 ['FtsZ division machinery','Cytoplasmic division site','Tubulin-related bacterial division protein; do not label it a eukaryotic microtubule. Atomic T. denticola model not assigned.','TD_MOTOR','internal']
];
export const extraSources=Object.values(extraProteins).flat().concat(tdMachinery).filter(p=>p[5]).map(p=>({id:p[3],url:'https://www.rcsb.org/structure/'+p[5],depiction:p[2],evidence:'Coordinate entry and construct identity; exact cell-level placement remains separate.'})).concat([
 {id:'EXT_9P47',url:'https://www.rcsb.org/structure/9P47',depiction:'T. denticola OppA peptide-binding complex; no surface anchor reconstructed.'},
 {id:'EXT_6KMF',url:'https://www.rcsb.org/structure/6KMF',depiction:'Polymerized P. gingivalis FimA segment; native attachment base absent.'},
 {id:'TD_CFP',url:'https://pmc.ncbi.nlm.nih.gov/articles/PMC94976/',depiction:'CfpA forms distinct internal cytoplasmic filaments; not flagellin or tubulin.'},
 {id:'TD_MOTOR',url:'https://pmc.ncbi.nlm.nih.gov/articles/PMC2519799/',depiction:'Tomography supports internal flagellar bundles, polar basal bodies and division-site organization.'},
 {id:'PM_TMP',url:'https://pmc.ncbi.nlm.nih.gov/articles/PMC13080454/',depiction:'Reported P. micra surface-associated TmpC; public accession verified separately, atomic coordinate download unavailable.'}
]);
export const placementIds=['2OMF','2FCP','4J7C','2B4L','5OF1','7AD3','7VH5'];
export const fitIds=['2OMF','2FCP','6SLJ','8BYS'];
export const assemblyIds=['9PIN','6KMF','8AUR','5W3H'];
