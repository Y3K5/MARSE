import * as T from 'three';
export function sampledLayer(Tdata,color,scale,alpha=.45){const {mid,normals}=Tdata,nu=64,nv=33,coords=new Float32Array(mid.flat().map(x=>x/scale)),indices=[];
 // Deliberate cutaway: omit a 90-degree sector; no atomistic surface claim.
 for(let u=0;u<nu;u++){if(u>=8&&u<24)continue;for(let v=0;v<nv-1;v++){const a=u*nv+v,b=((u+1)%nu)*nv+v;indices.push(a,b,a+1,b,b+1,a+1);}}
 const g=new T.BufferGeometry();g.setAttribute('position',new T.BufferAttribute(coords,3));g.setIndex(indices);g.computeVertexNormals();return new T.Mesh(g,new T.MeshStandardMaterial({color,side:T.DoubleSide,transparent:true,opacity:alpha,roughness:.43,metalness:.13,depthWrite:false}));}
export function machineryCell(t,detail,mesh,ball,tube,label,parent){
 const g=new T.Group();parent.add(g);
 if(t.id==='td'){
  const center=u=>[(u-.5)*8,Math.sin(u*Math.PI*5)*.55,Math.cos(u*Math.PI*5)*.55];
  const pts=Array.from({length:130},(_,i)=>new T.Vector3(...center(i/129)));
  const curve=new T.CatmullRomCurve3(pts);
  for(const [r,color,alpha]of [[.14,'#b6c78b',.25],[.105,'#b5a4ce',.45]]){const m=mesh(new T.TubeGeometry(curve,160,r,24,false),color,[0,0,0],alpha,g);}
  const f=Array.from({length:150},(_,i)=>{const u=i/149,p=center(u),a=u*30;return[p[0],p[1]+Math.cos(a)*.12,p[2]+Math.sin(a)*.12]});const filament=tube(f,.009,'#dfb389',g);filament.userData.component='Flagellar filament segment';
  const cfp=Array.from({length:110},(_,i)=>{const p=center(i/109);return[p[0],p[1]-.065,p[2]+.04]});const ribbon=tube(cfp,.006,'#8cbdcf',g);ribbon.userData.component='CfpA ribbon';
  for(const u of [.03,.97]){const p=center(u);const motor=new T.Group();motor.position.set(...p);g.add(motor);for(let j=0;j<3;j++){const tor=mesh(new T.TorusGeometry(.06+j*.008,.012,8,28),'#c4b3d5',[0,0,j*.03],1,motor);tor.userData.component='Flagellar motor + hook';}ball([p[0],p[1],p[2]+.08],.035,'#e6b58c',g);}
  // Division-site ring is a schematic contextual marker, not a polymer model.
  const p=center(.5);const ring=mesh(new T.TorusGeometry(.085,.009,8,40),'#d6b6b7',p,1,g);ring.rotation.y=Math.PI/2;ring.userData.component='FtsZ division machinery';
  label('Outer sheath · transparent cutaway',[0,1.3,0]);label('Periplasmic flagellar filament',[1.8,.9,0]);label('CfpA ribbon · inside inner membrane',[-1.8,-1.1,0]);label('Polar motor / hook · schematic',[-3.3,.6,0]);label('Division site · FtsZ context',[0,-1,0]);return g;
 }
 return null;
}
