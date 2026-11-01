import * as T from 'three';
import {recordStyle,edgePath} from './network-model.js';
export function makeNetworkEdge(record,layout,index,{emphasized=false,subdued=false,scale=1}={}){
 const style=recordStyle(record),points=edgePath(record,layout,index).map(p=>new T.Vector3(...p)),curve=new T.CatmullRomCurve3(points,false,'centripetal'),root=new T.Group(),across=layout[record.a].species!==layout[record.b].species,opacity=emphasized?1:subdued?.16:across?.96:.38,emphasisScale=emphasized?Math.max(1,scale*1.3):1;
 if(style.pattern==='solid')root.add(new T.Mesh(new T.TubeGeometry(curve,100,(across?.24:.10)*emphasisScale,across?10:6,false),new T.MeshStandardMaterial({color:style.color,roughness:.6,transparent:true,opacity})));
 else{const line=new T.Line(new T.BufferGeometry().setFromPoints(curve.getPoints(100)),new T.LineDashedMaterial({color:style.color,transparent:true,opacity,dashSize:style.pattern==='dot'?.35:1.8,gapSize:style.pattern==='dot'?.8:1.15}));line.computeLineDistances();root.add(line);}
 if(across&&style.pattern!=='solid'){
  const length=curve.getLength(),dash=style.pattern==='dot'?.65:2.5,gap=style.pattern==='dot'?1.05:1.6,count=Math.ceil(length/(dash+gap)),mesh=new T.InstancedMesh(new T.CylinderGeometry(.18,.18,1,8),new T.MeshStandardMaterial({color:style.color,roughness:.5,transparent:true,opacity}),count),o=new T.Object3D(),axis=new T.Vector3(0,1,0);
  for(let i=0;i<count;i++){const a=curve.getPointAt(Math.min(1,i*(dash+gap)/length)),b=curve.getPointAt(Math.min(1,(i*(dash+gap)+dash)/length)),direction=b.clone().sub(a);o.position.copy(a).add(b).multiplyScalar(.5);o.quaternion.setFromUnitVectors(axis,direction.clone().normalize());o.scale.set(1,Math.max(.0001,direction.length()),1);o.updateMatrix();mesh.setMatrixAt(i,o.matrix);}mesh.userData.networkDashSegments=true;root.add(mesh);
 }
 const marker=new T.Mesh(new T.SphereGeometry((across?.85:.48)*emphasisScale,16,12),new T.MeshStandardMaterial({color:style.color,roughness:.4,emissive:style.color,emissiveIntensity:.12}));marker.position.copy(curve.getPoint(.5));root.add(marker);
 const proxy=new T.Mesh(new T.TubeGeometry(curve,55,.72,5,false),new T.MeshBasicMaterial({transparent:true,opacity:0,depthWrite:false,color:style.color}));proxy.userData.pickProxy=true;root.add(proxy);
 if(record.direction==='directed'&&record.a!==record.b){const pos=curve.getPoint(.86),tangent=curve.getTangent(.86).normalize();let cap;
  if(style.key==='inhibition'){cap=new T.Mesh(new T.BoxGeometry(across?2:1.4,.23,.23),new T.MeshBasicMaterial({color:style.color}));cap.quaternion.setFromUnitVectors(new T.Vector3(0,1,0),tangent);}else{cap=new T.Mesh(new T.ConeGeometry(across?.58:.4,across?1.7:1.2,14),new T.MeshStandardMaterial({color:style.color}));cap.quaternion.setFromUnitVectors(new T.Vector3(0,1,0),tangent);}cap.position.copy(pos);root.add(cap);
 }
 root.userData.networkRecord={id:record.id,a:record.a,b:record.b,kind:style.key,level:record.level,relation:record.relation,direction:record.direction,start:[...layout[record.a].position],end:[...layout[record.b].position],selfLoop:record.a===record.b,across,supportedHighlight:emphasized};
 root.traverse(o=>{o.castShadow=false;o.receiveShadow=false;if(o.material&&!o.userData.pickProxy){o.userData.baseNetworkColor=style.color;o.userData.baseNetworkOpacity=o.material.opacity;}});return root;
}
export function highlightNetworkEdges(edges,selection){
 for(const edge of edges){const r=edge.userData.networkRecord,active=selection?.kind==='record'?r.id===selection.id:selection?.kind==='protein'||selection?.kind==='context'?r.a===selection.id||r.b===selection.id:false,hasSelection=['record','protein','context'].includes(selection?.kind);
  edge.traverse(o=>{if(!o.material||o.userData.pickProxy)return;const color=active?'#f1d49a':o.userData.baseNetworkColor;o.material.color?.set(color);o.material.transparent=true;o.material.opacity=active?1:hasSelection?.24:o.userData.baseNetworkOpacity;});
 }
}
