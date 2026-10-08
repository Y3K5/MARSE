import * as T from 'three';
import type {TissueId} from './anatomy';
// All patterns are deterministic artistic surface cues, not measured histology.
function texture(kind:'gum'|'bone'|'enamel'|'dentin'|'fiber'|'pulp',section=false){
 const size=1024,c=document.createElement('canvas');c.width=c.height=size;const ctx=c.getContext('2d')!,pixels=ctx.createImageData(size,size);
 let seed=813;const random=()=>{seed=(1664525*seed+1013904223)>>>0;return seed/4294967296;};
 for(let y=0;y<size;y++)for(let x=0;x<size;x++){
  const noise=random()-.5,fold=Math.sin(x*.017+Math.sin(y*.013)*2)+.38*Math.sin(y*.054+x*.024),i=(y*size+x)*4;
  let r=255,g=255,b=255;
  if(kind==='gum'){const v=fold*3+noise*9;r=176+v;g=88+v*.8;b=88+v*.7;}
  else if(kind==='fiber'){const v=3*Math.sin(x*.045+Math.sin(y*.012)*11)+noise*10;r=176+v;g=105+v;b=95+v;}
  else if(kind==='pulp'){const v=fold*9+noise*13;r=128+v;g=36+v*.5;b=48+v*.5;}
  else if(kind==='enamel'){const v=.65*Math.sin(x*.11+Math.sin(y*.007)*5)+noise*1.3;r=247+v;g=242+v;b=226+v;}
  else if(kind==='dentin'){const angle=Math.atan2(y-512,x-512),v=9*Math.sin(angle*275)+noise*9;r=218+v;g=183+v;b=123+v;}
  else {const v=fold*3+noise*7;r=157+v;g=143+v;b=119+v;}
  pixels.data[i]=r;pixels.data[i+1]=g;pixels.data[i+2]=b;pixels.data[i+3]=255;
 }
 ctx.putImageData(pixels,0,0);
 if(kind==='gum'){
  // Tiny stippling and faint capillary-like surface marks, purely schematic art.
  for(let i=0;i<6800;i++){const x=random()*size,y=random()*size,r=.4+random()*1.1;ctx.fillStyle='rgba(96,31,43,.14)';ctx.beginPath();ctx.ellipse(x,y,r,r*.65,0,0,Math.PI*2);ctx.fill();}
 }else if(kind==='bone'&&section){
  // Irregular connected struts rather than a repeated pore tile.
  ctx.fillStyle='#80694f';ctx.fillRect(0,0,size,size);
  const step=72,side=Math.ceil(size/step)+2,sites:Array<[number,number]>=[];
  for(let yy=-1;yy<side-1;yy++)for(let xx=-1;xx<side-1;xx++)sites.push([(xx+.06+random()*.88)*step,(yy+.06+random()*.88)*step]);
  ctx.lineCap='round';
  const strut=(a:[number,number],b:[number,number])=>{
   const dx=b[0]-a[0],dy=b[1]-a[1],bend=(random()-.5)*22,w=9+random()*9;
   ctx.beginPath();ctx.moveTo(...a);ctx.quadraticCurveTo((a[0]+b[0])/2-dy/step*bend,(a[1]+b[1])/2+dx/step*bend,...b);
   ctx.strokeStyle='#b39b75';ctx.lineWidth=w+5;ctx.stroke();ctx.strokeStyle='#d6c09a';ctx.lineWidth=w;ctx.stroke();
  };
  for(let y=0;y<side-1;y++)for(let x=0;x<side-1;x++){
   const p=sites[y*side+x];if(random()<.82)strut(p,sites[y*side+x+1]);if(random()<.86)strut(p,sites[(y+1)*side+x]);if(random()<.24)strut(p,sites[(y+1)*side+x+1]);
  }
 }
 const t=new T.CanvasTexture(c);t.colorSpace=T.SRGBColorSpace;t.wrapS=t.wrapT=T.RepeatWrapping;t.anisotropy=4;return t;
}
export function pocketMaterials():Record<TissueId,T.MeshPhysicalMaterial>{
 const gum=texture('gum'),bone=texture('bone'),enamel=texture('enamel'),dentin=texture('dentin'),pulp=texture('pulp');
 const flesh=(map:T.Texture,color='#ffffff')=>new T.MeshPhysicalMaterial({color,map,roughness:.48,clearcoat:.18,clearcoatRoughness:.34,bumpMap:map,bumpScale:.018,sheen:.12,sheenColor:'#d59080',sheenRoughness:.8,side:T.DoubleSide});
 return {
 enamel:new T.MeshPhysicalMaterial({color:'#ffffff',map:enamel,roughness:.27,clearcoat:.48,clearcoatRoughness:.18,ior:1.55,bumpMap:enamel,bumpScale:.003,side:T.DoubleSide}),
 dentin:new T.MeshPhysicalMaterial({color:'#ffffff',map:dentin,roughness:.68,side:T.DoubleSide}),
 pulp:flesh(pulp),cementum:new T.MeshPhysicalMaterial({color:'#bdab83',roughness:.8,side:T.DoubleSide}),
 pdl:new T.MeshPhysicalMaterial({color:'#a498b3',roughness:.72,side:T.DoubleSide}),gingiva:flesh(gum),
 bone:new T.MeshPhysicalMaterial({color:'#aa987d',map:bone,roughness:.85,bumpMap:bone,bumpScale:.045,side:T.DoubleSide}),
 vessels:flesh(pulp,'#b9423f'),nerve:new T.MeshPhysicalMaterial({color:'#caa259',roughness:.55,side:T.DoubleSide}),
 palate:flesh(gum),tongue:flesh(gum,'#d58887')
 };
}
export function boneSectionTexture(){return texture('bone',true);}
export function connectiveSectionTexture(){return texture('fiber',true);}
