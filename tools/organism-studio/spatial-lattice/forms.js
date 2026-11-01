import * as T from 'three';
export const mat=(color,options={})=>new T.MeshStandardMaterial({color,roughness:.53,metalness:.07,...options});
const sphere=new T.SphereGeometry(1,28,20);
export function ball(parent,r,pos,color,options={}){const m=new T.Mesh(sphere,mat(color,options));m.scale.setScalar(r);m.position.set(...pos);parent.add(m);return m;}
export function box(parent,size,pos,color,bevel=.1){const mesh=new T.Mesh(new T.BoxGeometry(...size,2,2,2),mat(color));mesh.position.set(...pos);parent.add(mesh);return mesh;}
export function tube(parent,points,r,color,segments=90,options={}){const curve=new T.CatmullRomCurve3(points.map(p=>new T.Vector3(...p)));const mesh=new T.Mesh(new T.TubeGeometry(curve,segments,r,8,false),mat(color,options));parent.add(mesh);return mesh;}
export function rod(parent,a,b,r,color,options={}){const av=new T.Vector3(...a),bv=new T.Vector3(...b),v=bv.clone().sub(av);const mesh=new T.Mesh(new T.CylinderGeometry(r,r,v.length(),20),mat(color,options));mesh.position.copy(av.add(bv).multiplyScalar(.5));mesh.quaternion.setFromUnitVectors(new T.Vector3(0,1,0),v.normalize());parent.add(mesh);return mesh;}
function instanceBalls(parent,positions,r,color){const mesh=new T.InstancedMesh(new T.SphereGeometry(1,10,8),mat(color),positions.length),temp=new T.Object3D();positions.forEach((p,i)=>{temp.position.set(...p);temp.scale.setScalar(typeof r==='function'?r(i):r);temp.updateMatrix();mesh.setMatrixAt(i,temp.matrix);});parent.add(mesh);return mesh;}
export function microbe(id,{interior=false,quality='high'}={}){
 const g=new T.Group(),inner=new T.Group();g.add(inner);const color={pg:'#d59181',tf:'#d5b170',td:'#77bdac'}[id];const opacity=interior?.24:1;
 if(id==='td'){
  const pts=[];for(let i=0;i<=240;i++){const t=i/240,rad=2.0*Math.sin(Math.PI*t)**.45;pts.push([(t-.5)*30,Math.sin(t*Math.PI*10)*rad,Math.cos(t*Math.PI*10)*rad]);}
  const helix=tube(g,pts,.58,color,quality==='high'?400:240,{transparent:interior,opacity,roughness:.4});helix.rotation.z=.35;
  if(interior){tube(inner,pts,.2,'#dcd7a5',240);inner.rotation.z=.35;}
 }else{
  const pg=id==='pg',length=pg?12:22,radius=pg?5:3.1;
  let geometry;
  if(pg){geometry=new T.CapsuleGeometry(radius,length-radius*2,20,quality==='high'?56:28);geometry.rotateZ(Math.PI/2);}else{
   const points=[];for(let i=0;i<=80;i++){const t=i/80;points.push(new T.Vector2(radius*Math.sin(Math.PI*t)**.65,(t-.5)*length));}geometry=new T.LatheGeometry(points,quality==='high'?64:32);geometry.rotateZ(Math.PI/2);
  }
  const body=new T.Mesh(geometry,mat(color,{transparent:interior,opacity,roughness:.6}));g.add(body);
  const positions=[],n=quality==='high'?340:150;
  for(let i=0;i<n;i++){const t=(i+.5)/n,x=(t-.5)*(pg?11:21),a=i*2.399963;let r=pg?Math.sqrt(Math.max(.2,radius*radius-Math.max(0,Math.abs(x)-1)**2)):radius*Math.sin(Math.PI*t)**.65;positions.push([x,Math.cos(a)*(r+.13),Math.sin(a)*(r+.13)]);}
  const coat=instanceBalls(g,positions,pg?.26:.22,pg?'#dea494':'#e0c489');coat.material.transparent=interior;coat.material.opacity=interior?.22:1;
  if(pg){for(let i=0;i<48;i++){const x=Math.sin(i*4.21)*4.4,a=i*2.39996,r=4.5,ext=1.8+(i%5)*.3;const p=[x,Math.cos(a)*r,Math.sin(a)*r];tube(g,[p,[x+.3,Math.cos(a)*(r+ext*.5),Math.sin(a)*(r+ext*.5)],[x+.8,Math.cos(a+.1)*(r+ext),Math.sin(a+.1)*(r+ext)]],.07,'#d8aa97',10,{transparent:interior,opacity:interior?.3:.85});}}
  if(interior){const dna=[];for(let i=0;i<120;i++){const a=i*.32;dna.push([Math.sin(a*.31)*length*.26,Math.sin(a)*1.5,Math.cos(a)*1.5]);}tube(inner,dna,.13,'#e1cf9f',180);instanceBalls(inner,Array.from({length:80},(_,i)=>[Math.sin(i*2.1)*length*.32,Math.cos(i*3.7)*radius*.65,Math.sin(i*4.2)*radius*.65]),.24,'#add5ce');}
 }
 g.rotation.set(id==='tf'?.25:.1,id==='pg'?-.4:.15,id==='tf'?-.5:.08);return g;
}
export function makePen({explode=0,color='#749cbd'}={}){
 const root=new T.Group(),parts=[];
 const part=(name,x,fn)=>{const g=new T.Group();g.position.x=x;fn(g);root.add(g);parts.push({name,object:g});};
 part('Barrel',0,g=>{const m=new T.Mesh(new T.CylinderGeometry(.52,.52,9,48,1,true),mat(color,{roughness:.32,metalness:.35}));m.rotation.z=Math.PI/2;g.add(m);for(let i=0;i<9;i++)rod(g,[-3+i*.18,-.2,-.44],[-3+i*.18,.2,-.44],.025,'#b7c6ce');});
 part('Ink refill',explode*2.5,g=>{g.position.y=explode*3;rod(g,[-5.7,0,0],[4.5,0,0],.12,'#bac5c8');rod(g,[-3.8,0,0],[3.8,0,0],.14,'#38505d');});
 part('Grip',-5-explode*2.2,g=>{const m=new T.Mesh(new T.CylinderGeometry(.43,.52,1.5,40),mat('#3b474f'));m.rotation.z=Math.PI/2;g.add(m);for(let i=0;i<7;i++){const r=new T.Mesh(new T.TorusGeometry(.49,.025,8,32),mat('#76878d'));r.rotation.y=Math.PI/2;r.position.x=(i-3)*.18;g.add(r);}});
 part('Tip',-6.2-explode*4.4,g=>{const m=new T.Mesh(new T.ConeGeometry(.41,1,40),mat('#c5c7c3',{metalness:.85,roughness:.23}));m.rotation.z=Math.PI/2;g.add(m);ball(g,.075,[-.55,0,0],'#d3d4c9',{metalness:.9,roughness:.18});});
 part('End cap',5+explode*4,g=>{rod(g,[-.5,0,0],[.5,0,0],.53,color,{metalness:.4,roughness:.3});box(g,[3,.1,.18],[-1,.6,0],'#bdc8ce');});
 return {root,parts};
}
export function desk(){
 const g=new T.Group(),top=box(g,[84,2,48],[0,-1,0],'#ffffff');for(const x of [-35,35])for(const z of [-18,18])box(g,[2.2,35,2.2],[x,-19,z],'#3f5059');
 const textureCanvas=document.createElement('canvas');textureCanvas.width=512;textureCanvas.height=128;const ctx=textureCanvas.getContext('2d'),image=ctx.createImageData(512,128);for(let y=0;y<128;y++)for(let x=0;x<512;x++){const i=(y*512+x)*4,v=Math.sin(y*.8+Math.sin(x*.012)*1.4)*3+Math.sin(y*.23+x*.004)*2;image.data.set([151+v,115+v,78+v,255],i);}ctx.putImageData(image,0,0);const texture=new T.CanvasTexture(textureCanvas);texture.colorSpace=T.SRGBColorSpace;top.material.map=texture;top.material.roughness=.8;top.material.needsUpdate=true;
 const notebook=new T.Group();notebook.position.set(-18,0,-2);box(notebook,[19,1,25],[0,.55,0],'#bbd0c1');box(notebook,[18.5,.7,24.5],[0,1.1,0],'#e5dfcc');for(let i=0;i<9;i++)box(notebook,[15,.025,.025],[0,1.47,-9+i*2],'#b9bbae');g.add(notebook);
 const mug=new T.Group(),outer=new T.Mesh(new T.CylinderGeometry(3.5,3,7,48,1,true),mat('#bac3b3',{side:T.DoubleSide,roughness:.4}));outer.position.y=3.6;mug.add(outer);const rim=new T.Mesh(new T.TorusGeometry(3.5,.17,12,64),mat('#d3d9cb'));rim.rotation.x=Math.PI/2;rim.position.y=7.1;mug.add(rim);const drink=new T.Mesh(new T.CircleGeometry(3.3,48),mat('#382c25',{roughness:.2}));drink.rotation.x=-Math.PI/2;drink.position.y=6.4;mug.add(drink);const handle=new T.Mesh(new T.TorusGeometry(2,.45,16,40,Math.PI*1.7),mat('#bac3b3'));handle.position.set(3.5,3.6,0);mug.add(handle);mug.position.set(24,0,-10);g.add(mug);
 const pen=makePen();pen.root.position.set(4,.65,7);pen.root.rotation.y=-.3;g.add(pen.root);const pen2=makePen({color:'#d1ab78'});pen2.root.position.set(10,.65,10);pen2.root.rotation.y=.4;g.add(pen2.root);
 return {root:g,pen:pen.root,pen2:pen2.root,mug,notebook};
}
export function cell({cutaway=false}={}){
 const root=new T.Group(),shell=ball(root,7,[0,0,0],'#86bfb3',{transparent:true,opacity:cutaway?.14:.62,roughness:.38,depthWrite:false});shell.scale.set(7,5.04,6.3);
 ball(root,2.2,[-1,0,0],'#cfadba',{roughness:.55});ball(root,.65,[-1.5,.5,.6],'#e6c8b7');
 for(let i=0;i<7;i++){const a=i*2.4,x=Math.cos(a)*4,y=Math.sin(a*1.3)*2,z=Math.sin(a)*4;const mito=ball(root,1,[x,y,z],'#deb580');mito.scale.set(1.5,.6,.7);mito.rotation.y=a;const pts=[];for(let j=0;j<30;j++)pts.push([x+(j/29-.5)*2.3,y+Math.sin(j*.8)*.24,z+.56]);tube(root,pts,.06,'#9f795b',45);}
 for(let i=0;i<6;i++){const points=[];for(let j=0;j<40;j++){const a=j/39*Math.PI*1.5;points.push([Math.cos(a)*(2.8+i*.22)-1,Math.sin(a)*1.5,1.1+i*.24]);}tube(root,points,.12,'#c3b1c7',60);}
 instanceBalls(root,Array.from({length:120},(_,i)=>[Math.sin(i*2.43)*5,Math.sin(i*4.21)*2.4,Math.cos(i*2.9)*4]),.11,'#d7d2b4');return root;
}
export function tissue(){const root=new T.Group(),cells=[];for(let x=-2;x<=2;x++)for(let z=-1;z<=1;z++){const c=cell();c.scale.set(.88,.7,.88);c.position.set(x*12.5,(Math.sin(x*2+z)*.8),z*11);root.add(c);cells.push(c);}box(root,[70,.6,40],[0,-4.8,0],'#6a7785');return {root,cells};}
