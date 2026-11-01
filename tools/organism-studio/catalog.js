// Each source ID resolves in sources.json; geometry remains a separately labeled schematic.
export const proteins={
 td:[
 ['Msp / MOSP','Outer sheath and periplasm','Some regions exposed; stable periplasmic forms also occur.','TD2,TD3','surface'],
 ['Dentilisin','Outer-sheath complex','PrtP cleaves proteins; PrcA1/PrcB are lipid-anchored partners. PrcA2 is less accessible on intact wild-type cells.','TD2','surface'],
 ['FhbB','Surface lipoprotein','Binds host factor H. The experimental structure is a recombinant derivative, not a native anchored reconstruction.','TD4','surface','3QZ0'],
 ['OppA','Reported surface-associated protein','Host-protein binding reported. Exact placement and anchoring are not reconstructed.','TD5','uncertain'],
 ['FlaA + FlaB1–3','Periplasmic flagella','FlaA coats the filament; FlaB proteins form its core. Beneath the outer membrane, not exposed on an intact cell.','TD1','internal']],
 pg:[
 ['FimA fimbriae','External polymer','Major attachment fiber; FimC/D/E accessory components; FimB at the base.','PG1,PG3','surface'],
 ['Mfa1 fimbriae','External polymer','Minor attachment fiber; Mfa3/4/5 accessory components; Mfa2 basal component.','PG2,PG3','surface'],
 ['RgpA / RgpB / Kgp','Surface-associated proteases','RgpA/B cleave after arginine, Kgp after lysine. Extracellular and vesicle-associated forms occur.','PG4','surface'],
 ['RagAB','Outer-membrane transporter','RagA forms a membrane-spanning transporter; RagB forms an external lid. Peptide acquisition.','PG5','surface','6SLJ'],
 ['HmuR + HmuY','Heme-acquisition system','HmuR is an outer-membrane receptor; HmuY binds heme outside. No exact coupled placement inferred.','PG6','surface'],
 ['Pgm6 + Pgm7','OmpA-like complex','Envelope organization; parts face outside and parts lie in the periplasm. Cartoon does not specify native stoichiometry.','PG7','surface']],
 tf:[
 ['TfsA + TfsB','S-layer glycoproteins','Form the repeating outer coat. Glycans and lattice contacts can mask protein surfaces.','TF1','surface'],
 ['BspA','Surface-associated adhesin','Supports attachment and bacterial interactions; position relative to the S-layer unresolved here.','TF2','uncertain'],
 ['NanH','Surface-associated and secreted enzyme','Removes terminal sialic acids. Experimental coordinates do not establish a native membrane anchor.','TF3','surface','7QZ3'],
 ['NanO + NanU','Proposed nutrient-acquisition system','NanO is a proposed outer-membrane transporter; NanU binds sialic acid. Intact complex and S-layer access not reconstructed.','TF4,TF5','uncertain']]
};
export function validateFrames(d){
 if(!d || !Number.isInteger(d.width)||!Number.isInteger(d.height)||d.width<1||d.height<1||d.width*d.height>100000)throw Error('Expected a MARSE width/height grid of 1–100,000 cells.');
 if(!Array.isArray(d.species)||!d.species.length||d.species.length>50||new Set(d.species).size!==d.species.length||d.species.some(s=>typeof s!=='string'||!s.length||s.length>100||['__proto__','constructor','prototype'].includes(s)))throw Error('Expected unique named species.');
 if(!Array.isArray(d.frames)||!d.frames.length||d.frames.length>1000)throw Error('Expected 1–1,000 saved frames.');
 if(d.width*d.height*d.species.length*d.frames.length>10000000)throw Error('Frame collection exceeds the local display limit.');
 for(const f of d.frames){if(!Number.isFinite(f.time_h))throw Error('Every frame needs a finite time_h.');for(const s of d.species){const a=f.species?.[s];if(!Array.isArray(a)||a.length!==d.height||a.some(r=>!Array.isArray(r)||r.length!==d.width||r.some(v=>!Number.isFinite(v)||v<0)))throw Error('Invalid or negative species density grid.');}}
 return d;
}
export function particles(c,taxa,max=240){
 const out=[]; const totals=c.biomass.map(a=>a.reduce((s,x)=>s+x,0));const total=totals.reduce((a,b)=>a+b,0);
 if(!total)return out;
 for(let k=0;k<taxa.length;k++){const n=Math.round(max*totals[k]/total);if(!totals[k])continue;let cursor=0,sum=c.biomass[k][0];for(let j=0;j<n;j++){const goal=totals[k]*(j+.5)/n;while(sum<goal&&cursor<c.biomass[k].length-1)sum+=c.biomass[k][++cursor];out.push({taxon:taxa[k].id,x:(cursor%c.nx+.5)/c.nx,y:(Math.floor(cursor/c.nx)+.5)/c.ny,z:Math.sin(j*12.7+k)*.5});}}
 return out;
}
