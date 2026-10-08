import * as T from 'three';
import type {TissueId} from './anatomy';
// All patterns are deterministic artistic surface cues, not measured histology.
function texture(kind:'gum'|'bone'|'enamel'|'dentin'|'fiber'|'pulp',section=false){
 const size=1024,c=document.createElement('canvas');c.width=c.height=size;const ctx=c.getContext('2d')!,pixels=ctx.createImageData(size,size);
 let seed=813;const random=()=>{seed=(1664525*seed+1013904223)>>>0;return seed/4294967296;};
 for(let y=0;y<size;y++)for(let x=0;x<size;x++){
  const noise=random()-.5,fold=Math.sin(x*.017+Math.sin(y*.013)*2)+.38*Math.sin(y*.054+x*.024),i=(y*size+x)*4;
  let r=255,g=255,b=255;
  // Natural coral soft tissue, warm ivory-tan dentin and pale cortical bone.
  if(kind==='gum'){const v=fold*2.5+noise*7;r=216+v;g=122+v*.8;b=114+v*.7;}
  else if(kind==='fiber'){const v=4*Math.sin(x*.045+Math.sin(y*.012)*11)+noise*8;r=236+v;g=192+v;b=178+v;}
  else if(kind==='pulp'){const v=fold*6+noise*9;r=196+v;g=104+v*.6;b=100+v*.6;}
  else if(kind==='enamel'){const v=.55*Math.sin(x*.11+Math.sin(y*.007)*5)+noise*1.1;r=250+v;g=245+v;b=232+v;}
  else if(kind==='dentin'){const angle=Math.atan2(y-512,x-512),v=5*Math.sin(angle*275)+noise*6;r=232+v;g=206+v;b=166+v;}
  else {const v=fold*2+noise*6;r=228+v;g=222+v;b=208+v;}
  pixels.data[i]=r;pixels.data[i+1]=g;pixels.data[i+2]=b;pixels.data[i+3]=255;
 }
 ctx.putImageData(pixels,0,0);
 if(kind==='gum'){
  // Tiny stippling and faint capillary-like surface marks, purely schematic art.
  for(let i=0;i<6800;i++){const x=random()*size,y=random()*size,r=.4+random()*1.1;ctx.fillStyle='rgba(150,62,62,.10)';ctx.beginPath();ctx.ellipse(x,y,r,r*.65,0,0,Math.PI*2);ctx.fill();}
 }else if(kind==='bone'&&!section){
  // Faint cortical porosity and vascular foramina.
  for(let i=0;i<1400;i++){const x=random()*size,y=random()*size,r=.5+random()*1.6;ctx.fillStyle=`rgba(120,96,70,${.08+random()*.12})`;ctx.beginPath();ctx.ellipse(x,y,r,r*(.5+random()*.5),random()*3,0,Math.PI*2);ctx.fill();}
 }else if(kind==='bone'&&section){
  // Irregular connected struts rather than a repeated pore tile.
  // Marrow spaces between finer, irregular ivory trabeculae.
  ctx.fillStyle='#8f6253';ctx.fillRect(0,0,size,size);
  for(let i=0;i<900;i++){const x=random()*size,y=random()*size,r=4+random()*14;ctx.fillStyle=`rgba(${120+random()*30},${70+random()*20},${58+random()*14},.35)`;ctx.beginPath();ctx.ellipse(x,y,r,r*(.6+random()*.4),random()*3,0,Math.PI*2);ctx.fill();}
  const step=46,side=Math.ceil(size/step)+2,sites:Array<[number,number]>=[];
  for(let yy=-1;yy<side-1;yy++)for(let xx=-1;xx<side-1;xx++)sites.push([(xx+.06+random()*.88)*step,(yy+.06+random()*.88)*step]);
  ctx.lineCap='round';
  const strut=(a:[number,number],b:[number,number])=>{
   const dx=b[0]-a[0],dy=b[1]-a[1],bend=(random()-.5)*20,w=4+random()*7;
   ctx.beginPath();ctx.moveTo(...a);ctx.quadraticCurveTo((a[0]+b[0])/2-dy/step*bend,(a[1]+b[1])/2+dx/step*bend,...b);
   ctx.strokeStyle='#c8b593';ctx.lineWidth=w+3;ctx.stroke();ctx.strokeStyle='#e6d9bd';ctx.lineWidth=w;ctx.stroke();
  };
  for(let y=0;y<side-1;y++)for(let x=0;x<side-1;x++){
   const p=sites[y*side+x];if(random()<.72)strut(p,sites[y*side+x+1]);if(random()<.76)strut(p,sites[(y+1)*side+x]);if(random()<.22)strut(p,sites[(y+1)*side+x+1]);
  }
 }
 const t=new T.CanvasTexture(c);t.colorSpace=T.SRGBColorSpace;t.wrapS=t.wrapT=T.RepeatWrapping;t.anisotropy=4;return t;
}
export function pocketMaterials():Record<TissueId,T.MeshPhysicalMaterial>{
 const gum=texture('gum'),bone=texture('bone'),enamel=texture('enamel'),dentin=texture('dentin'),pulp=texture('pulp');
 const flesh=(map:T.Texture,color='#ffffff')=>new T.MeshPhysicalMaterial({color,map,roughness:.5,clearcoat:.22,clearcoatRoughness:.38,bumpMap:map,bumpScale:.014,sheen:.25,sheenColor:'#f0a89a',sheenRoughness:.7,side:T.DoubleSide});
 return {
 // Ivory enamel with a soft, controlled sheen rather than a mirror-like clearcoat.
 enamel:new T.MeshPhysicalMaterial({color:'#ffffff',map:enamel,roughness:.36,clearcoat:.22,clearcoatRoughness:.32,specularIntensity:.55,ior:1.55,sheen:.18,sheenColor:'#fff4dc',sheenRoughness:.6,bumpMap:enamel,bumpScale:.003,side:T.DoubleSide}),
 dentin:new T.MeshPhysicalMaterial({color:'#ffffff',map:dentin,roughness:.7,side:T.DoubleSide}),
 pulp:flesh(pulp),cementum:new T.MeshPhysicalMaterial({color:'#cdb58c',roughness:.82,side:T.DoubleSide}),
 pdl:new T.MeshPhysicalMaterial({color:'#b08a96',roughness:.74,side:T.DoubleSide}),gingiva:flesh(gum),
 bone:new T.MeshPhysicalMaterial({color:'#d6cdb9',map:bone,roughness:.88,bumpMap:bone,bumpScale:.03,side:T.DoubleSide}),
 vessels:flesh(pulp,'#b9423f'),nerve:new T.MeshPhysicalMaterial({color:'#caa259',roughness:.55,side:T.DoubleSide}),
 palate:flesh(gum),tongue:flesh(gum,'#d58887')
 };
}
export function boneSectionTexture(){return texture('bone',true);}
export function connectiveSectionTexture(){return texture('fiber',true);}
