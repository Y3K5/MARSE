import * as T from 'three';
import {RoomEnvironment} from 'three/addons/environments/RoomEnvironment.js';
import type {TissueId} from './anatomy';
export function environment(renderer:T.WebGLRenderer):T.Texture {
 const pmrem=new T.PMREMGenerator(renderer),room=new RoomEnvironment();const map=pmrem.fromScene(room,.04).texture;room.dispose();pmrem.dispose();return map;
}
function texture(kind:'gum'|'bone'|'enamel'|'tongue'):T.CanvasTexture {
 const c=document.createElement('canvas');c.width=c.height=512;const ctx=c.getContext('2d')!;const data=ctx.createImageData(512,512);
 let seed=173;const random=()=>{seed=(1664525*seed+1013904223)>>>0;return seed/4294967296;};
 for(let y=0;y<512;y++)for(let x=0;x<512;x++){
  const n=random(),wave=Math.sin(x*(kind==='enamel'?.8:.11)+Math.sin(y*.07)*2),i=(y*512+x)*4;
  const v=kind==='enamel'?155+wave*14+n*18:kind==='bone'?185+n*40:125+n*55;
  data.data[i]=data.data[i+1]=data.data[i+2]=v;data.data[i+3]=255;
 }
 ctx.putImageData(data,0,0);
 if(kind!=='enamel')for(let i=0;i<(kind==='bone'?180:2100);i++){
  const x=random()*512,y=random()*512,r=kind==='bone'?1+random()*7:.3+random()*1.1;
  ctx.fillStyle=kind==='bone'?'#3a3a3a':'#bbbbbb';ctx.beginPath();ctx.ellipse(x,y,r,r*(.6+random()*.4),random()*3,0,Math.PI*2);ctx.fill();
 }
 const map=new T.CanvasTexture(c);map.wrapS=map.wrapT=T.RepeatWrapping;map.repeat.set(kind==='enamel'?1:3,3);return map;
}
export function materials():Record<TissueId,T.MeshPhysicalMaterial> {
 const gum=texture('gum'),bone=texture('bone'),enamel=texture('enamel'),tongue=texture('tongue');
 const flesh=(color:string)=>new T.MeshPhysicalMaterial({color,roughness:.36,clearcoat:.45,clearcoatRoughness:.18,bumpMap:gum,bumpScale:.045,side:T.DoubleSide});
 return {
 enamel:new T.MeshPhysicalMaterial({color:'#eee5ce',roughness:.2,metalness:0,clearcoat:1,clearcoatRoughness:.09,ior:1.55,bumpMap:enamel,bumpScale:.017,side:T.DoubleSide}),
 dentin:new T.MeshPhysicalMaterial({color:'#d4af72',roughness:.6,side:T.DoubleSide}),
 pulp:flesh('#ad4855'),cementum:new T.MeshPhysicalMaterial({color:'#cbb890',roughness:.7,side:T.DoubleSide}),
 pdl:flesh('#bc9785'),gingiva:flesh('#b76872'),
 bone:new T.MeshPhysicalMaterial({color:'#d2bf98',roughness:.9,map:bone,bumpMap:bone,bumpScale:.13,side:T.DoubleSide}),
 vessels:flesh('#b92f49'),nerve:new T.MeshPhysicalMaterial({color:'#dfc174',roughness:.45,side:T.DoubleSide}),
 palate:new T.MeshPhysicalMaterial({color:'#924858',roughness:.55,clearcoat:.18,clearcoatRoughness:.3,bumpMap:gum,bumpScale:.07,side:T.DoubleSide}),tongue:new T.MeshPhysicalMaterial({color:'#a94f60',roughness:.42,clearcoat:.42,clearcoatRoughness:.22,bumpMap:tongue,bumpScale:.10,side:T.DoubleSide})
 };
}
