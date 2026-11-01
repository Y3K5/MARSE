import * as T from './vendor/three.module.js';
// This is a teaching/design lattice of glyphs. It never invents atomic packing.
export function latticeDesign({protofilaments=12,rows=4,radiusNm=12}={}){
 if(!Number.isInteger(protofilaments)||![12,13].includes(protofilaments)||!Number.isInteger(rows)||rows<2||rows>8||!Number.isFinite(radiusNm)||radiusNm<8||radiusNm>20)throw Error('Invalid lattice design dimensions');
 const spacingNm=8,instances=[];
 for(let p=0;p<protofilaments;p++)for(let r=0;r<rows;r++){const angle=2*Math.PI*p/protofilaments;instances.push({id:`PF${p+1}-D${r+1}`,source:'5W3H',protofilament:p+1,row:r+1,center_nm:[radiusNm*Math.cos(angle),radiusNm*Math.sin(angle),(r-(rows-1)/2)*spacingNm],axis:[0,0,1]});}
 return{schema:'marse.tubulin_layout/1',protofilaments,rows,radiusNm,spacingNm,instances,units:'nm',source:'https://www.rcsb.org/structure/5W3H',architecture_source:'https://pmc.ncbi.nlm.nih.gov/articles/PMC5584162/',state:'user_design_glyph_layout',limits:['12/13 protofilament architectures are discussed in the source study; this displayed packing is not its reconstruction','Radius and 8 nm repeat are design presets, not fitted values','Glyphs link to the experimental tubulin dimer; atomic coordinates are not duplicated or physically packed']};
}
export function drawLattice(design,group,mesh,label,picks){
 const scale=5,shape=new T.SphereGeometry(1,20,14);
 for(const unit of design.instances){for(const [part,col]of [[-1,'#aacfc7'],[1,'#d5af8b']]){const p=unit.center_nm;const m=mesh(shape.clone(),col,[p[0]/scale,p[1]/scale,(p[2]+part*2)/scale],1,group);m.scale.set(.36,.36,.38);m.userData={sourceProtein:'5W3H',layoutUnit:unit.id};picks.push(m);}}
 label('α / β subunits · schematic glyphs',[0,0,design.rows*.85]);label('Select a unit to inspect its real fold',[0,0,-design.rows*.85]);
}
