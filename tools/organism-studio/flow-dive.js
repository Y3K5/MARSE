import * as T from 'three';

const palette = ['#f2bd8d','#8ccfd2','#bbabe9','#dfb1cc','#bad38b','#86b5e7'];
const up = new T.Vector3(0,1,0);
const material = color => new T.MeshStandardMaterial({color,roughness:.48,metalness:.04});
function dispose(node){node.traverse(o=>{o.geometry?.dispose();for(const m of [].concat(o.material||[]))m.dispose();});node.removeFromParent();}
function instances(geometry, color, points, parent, radius,weights){
 const mesh=new T.InstancedMesh(geometry,material(color),points.length),dummy=new T.Object3D();
 points.forEach((p,i)=>{dummy.position.copy(p);dummy.scale.setScalar(radius*(weights?.[i]??1));dummy.updateMatrix();mesh.setMatrixAt(i,dummy.matrix);});
 parent.add(mesh);return mesh;
}
function cylinders(pairs,radius,color,parent){
 const mesh=new T.InstancedMesh(new T.CylinderGeometry(1,1,1,8),material(color),pairs.length),dummy=new T.Object3D();
 pairs.forEach(([a,b],i)=>{const v=b.clone().sub(a);dummy.position.copy(a).add(b).multiplyScalar(.5);dummy.quaternion.setFromUnitVectors(up,v.clone().normalize());dummy.scale.set(radius,v.length(),radius);dummy.updateMatrix();mesh.setMatrixAt(i,dummy.matrix);});parent.add(mesh);return mesh;
}

// Connectivity is built from the displayed triangles, welded only at equal positions.
// No through-body ray samples: geodesic growth cannot jump to another coil or shell.
export function surfaceGraph(surface){
 const geometry=surface.geometry,p=geometry.attributes.position,n=geometry.attributes.normal;
 const normalMatrix=new T.Matrix3().getNormalMatrix(surface.matrixWorld),vertices=[],map=new Map(),ids=[];
 for(let i=0;i<p.count;i++){
  const point=new T.Vector3().fromBufferAttribute(p,i).applyMatrix4(surface.matrixWorld),key=point.toArray().map(v=>v.toFixed(6)).join(',');
  let id=map.get(key);if(id===undefined){id=vertices.length;map.set(key,id);vertices.push({p:point,n:new T.Vector3(),edges:new Map(),used:false});}
  ids[i]=id;if(n)vertices[id].n.add(new T.Vector3().fromBufferAttribute(n,i).applyMatrix3(normalMatrix).normalize());
 }
 const index=geometry.index,triangles=[];
 for(let i=0;i<(index?.count||p.count);i+=3){const tri=[0,1,2].map(k=>ids[index?index.getX(i+k):i+k]);if(new Set(tri).size<3)continue;triangles.push(tri);for(let k=0;k<3;k++){const a=vertices[tri[k]],b=vertices[tri[(k+1)%3]];a.used=b.used=true;const d=a.p.distanceTo(b.p);a.edges.set(tri[(k+1)%3],d);b.edges.set(tri[k],d);}}
 vertices.forEach(v=>v.n.normalize());return {vertices,triangles};
}

export class FlowDive{
 constructor({group,camera,controls,record,onStatus,organism,component,atoms=[]}){
  Object.assign(this,{group,camera,controls,record,onStatus,organism,component,rawAtoms:atoms});
  this.active=true;this.level=0;this.targetLevel=0;this.transition=null;
  this.start=camera.position.clone();this.startTarget=controls.target.clone();this.sourcePoints=[];this.displayPoints=[];this.representation='backbone';
  group.updateMatrixWorld(true);group.traverse(o=>{if(!this.surface&&o.isMesh&&o.geometry.attributes.normal&&(o.userData.cellBoundary||!organism&&!o.userData.component))this.surface=o;});
  if(!this.surface)throw Error('No inspectable body surface for this design.');
  this.bodySurfaces=[];group.traverse(o=>{if(o.isMesh&&o.geometry.attributes.normal&&(o.userData.cellBoundary||!organism&&!o.userData.component))this.bodySurfaces.push(o);});
  this.contextMeshes=[];group.traverse(o=>{if(o.isMesh&&!this.bodySurfaces.includes(o))this.contextMeshes.push(o);});
  this.graph=surfaceGraph(this.surface);
  const box=new T.Box3().setFromObject(this.surface),center=box.getCenter(new T.Vector3()),view=camera.position.clone().sub(center).normalize();
  let best=-Infinity,chosen=0;this.graph.vertices.forEach((v,i)=>{if(!v.used||v.n.lengthSq()<.5)return;const delta=v.p.clone().sub(center),score=delta.dot(view)-Math.abs(delta.x)*2-Math.abs(delta.y)*.2;if(score>best){best=score;chosen=i;}});
  this.root=new T.Group();this.root.name='selected-surface-region';this.root.userData={organism,component,source:record?.id,sha256:record?.sha256,placement:'linked inspection only; no native pose'};group.add(this.root);
  this.region=new T.Group();this.root.add(this.region);this.protein=new T.Group();this.protein.name='source-coordinate-inspection';this.root.add(this.protein);
  this.regionSerial=0;this.selectRegion({point:this.graph.vertices[chosen].p,vertex:chosen},false);
  this.buildProtein();
  // Whole-cell framing includes the linked cradle, while keeping the cell centre as target.
  const wholeBox=new T.Box3().setFromObject(group),direction=this.start.clone().sub(this.startTarget).normalize(),right=new T.Vector3().crossVectors(up,direction).normalize(),vertical=new T.Vector3().crossVectors(direction,right).normalize(),tv=Math.tan(camera.fov*Math.PI/360),th=tv*camera.aspect;let distance=1;
  for(const x of [wholeBox.min.x,wholeBox.max.x])for(const y of [wholeBox.min.y,wholeBox.max.y])for(const z of [wholeBox.min.z,wholeBox.max.z]){const p=new T.Vector3(x,y,z).sub(this.startTarget);distance=Math.max(distance,p.dot(direction)+Math.abs(p.dot(right))/th,p.dot(direction)+Math.abs(p.dot(vertical))/tv);}this.start.copy(this.startTarget).addScaledVector(direction,distance*1.1);
  controls.minDistance=this.patchRadius*.15;controls.maxDistance=Math.max(150,this.start.length()*2);controls.enabled=true;
 }
 selectRegion(hit,move=true){
  if(hit.object&&hit.object!==this.surface){this.surface=hit.object;this.graph=surfaceGraph(this.surface);}
  const {vertices,triangles}=this.graph;let seed=hit.vertex;
  if(seed===undefined){let d=Infinity;vertices.forEach((v,i)=>{const x=v.p.distanceToSquared(hit.point);if(v.used&&x<d){d=x;seed=i;}});}
  const vertex=vertices[seed];this.anchor=vertex.p.clone();this.normal=vertex.n.clone();
  this.root.position.copy(this.group.worldToLocal(this.anchor.clone()));
  const parentRotation=this.group.getWorldQuaternion(new T.Quaternion());this.root.quaternion.copy(parentRotation.invert()).multiply(new T.Quaternion().setFromUnitVectors(new T.Vector3(0,0,1),this.normal));this.root.updateMatrixWorld(true);
  const radii=[];for(const [id,d] of vertex.edges){const angle=vertex.n.angleTo(vertices[id].n);if(angle>.05)radii.push(d/(2*Math.sin(angle/2)));}radii.sort((a,b)=>a-b);
  const scale=this.surface.getWorldScale(new T.Vector3()),declared=this.surface.geometry.parameters?.radius;
  this.localRadius=declared?declared*Math.min(scale.x,scale.y,scale.z):(radii[Math.floor(radii.length*.25)]||new T.Box3().setFromObject(this.surface).getSize(new T.Vector3()).length()*.07);
  this.patchRadius=this.localRadius*.78;
  const distances=vertices.map(()=>Infinity);distances[seed]=0;const open=new Set([seed]);
  while(open.size){let id=-1,best=Infinity;for(const i of open)if(distances[i]<best){best=distances[i];id=i;}open.delete(id);if(best>this.patchRadius*2)continue;for(const [j,weight]of vertices[id].edges){const d=best+weight;if(d<distances[j]&&d<this.patchRadius*2.8){distances[j]=d;open.add(j);}}}
  while(this.region.children.length)dispose(this.region.children[0]);
  this.patchSamples=[];this.patchTriangles=[];const seen=new Set(),spacing=this.patchRadius/8;
  const inverse=this.root.getWorldQuaternion(new T.Quaternion()).invert();
  for(const tri of triangles){if(Math.min(...tri.map(i=>distances[i]))>this.patchRadius*2)continue;const vs=tri.map(i=>vertices[i]);
   const subdivisions=Math.min(32,Math.max(2,Math.ceil(Math.max(vs[0].p.distanceTo(vs[1].p),vs[1].p.distanceTo(vs[2].p),vs[2].p.distanceTo(vs[0].p))/spacing)));
   for(let i=0;i<=subdivisions;i++)for(let j=0;j<=subdivisions-i;j++){
    const weights=[i/subdivisions,j/subdivisions,1-(i+j)/subdivisions],distance=weights.reduce((s,w,k)=>s+(w?distances[tri[k]]*w:0),0);if(distance>this.patchRadius*2)continue;
    const world=new T.Vector3(),n=new T.Vector3();weights.forEach((w,k)=>{world.addScaledVector(vs[k].p,w);n.addScaledVector(vs[k].n,w);});n.normalize();const footprintDistance=world.distanceTo(this.anchor);if(footprintDistance>this.patchRadius)continue;
    if(n.dot(this.normal)<.25)continue;const key=world.toArray().map(v=>v.toFixed(5)).join(',');if(seen.has(key))continue;seen.add(key);
    this.patchSamples.push({p:this.root.worldToLocal(world.clone()),n:n.applyQuaternion(inverse),edge:T.MathUtils.smoothstep(1-footprintDistance/this.patchRadius,0,.25),triangle:tri,weights});
   }
  }
  const heads=[],inner=[],tails=[],weights=[];for(const {p,n,edge}of this.patchSamples){if(edge<.12)continue;const a=p.clone().addScaledVector(n,spacing*.13),b=p.clone().addScaledVector(n,-this.patchRadius*.105);heads.push(a);inner.push(b);tails.push([a,b]);weights.push(Math.max(.08,edge));}
  instances(new T.SphereGeometry(1,8,6),'#cee4b5',heads,this.region,spacing*.26,weights);
  instances(new T.SphereGeometry(1,8,6),'#8cbbbd',inner,this.region,spacing*.26,weights);
  cylinders(tails,spacing*.09,'#a4b59c',this.region);
  this.marker=new T.Mesh(new T.TorusGeometry(this.patchRadius*.14,spacing*.12,8,32),material('#f5d6a4'));this.marker.position.z=spacing*.2;this.region.add(this.marker);
  // The cradle is offset beyond the selected surface, never posed as an embedded component.
  this.protein.position.set(this.patchRadius*1.35,0,this.patchRadius*.75);
  const linkPoints=[new T.Vector3(0,0,spacing),new T.Vector3(this.patchRadius*.5,0,this.patchRadius*.75),this.protein.position.clone()];
  const link=new T.Line(new T.BufferGeometry().setFromPoints(linkPoints),new T.LineDashedMaterial({color:'#cba774',dashSize:this.patchRadius*.07,gapSize:this.patchRadius*.04}));link.computeLineDistances();this.region.add(link);
  this.regionSerial++;if(this.record&&this.factor){this.protein.scale.setScalar(this.factor);this.updateDisplayPoints();}if(move)this.setLevel(Math.max(1,this.targetLevel));
 }
 buildProtein(){
  if(!this.record)return;
  const rows=Object.entries(this.record.chains).flatMap(([chain,rs])=>rs.map(r=>({...r,chain})));this.residues=rows;
  const box=new T.Box3().setFromPoints(rows.map(r=>new T.Vector3(...r.p)));this.center=box.getCenter(new T.Vector3());this.sourceExtent=Math.max(...box.getSize(new T.Vector3()).toArray());
  // Fixed schematic magnification for this selection. It does not change with camera depth.
  this.factor=this.patchRadius*1.35/this.sourceExtent;this.protein.scale.setScalar(this.factor);
  const extents=box.getSize(new T.Vector3()).toArray(),axes=[0,1,2].sort((a,b)=>extents[b]-extents[a]),basis=[new T.Vector3(),new T.Vector3(),new T.Vector3()];basis.forEach((v,i)=>v.setComponent(axes[i],1));if(new T.Vector3().crossVectors(basis[0],basis[1]).dot(basis[2])<0)basis[2].negate();const rotation=new T.Matrix4().makeBasis(...basis).transpose();this.protein.quaternion.setFromRotationMatrix(rotation);
  this.backbone=new T.Group();this.protein.add(this.backbone);this.pickMeshes=[];
  let offset=0;for(const [i,[chain,rs]]of Object.entries(this.record.chains).entries()){
   const points=rs.map(r=>new T.Vector3(...r.p).sub(this.center)),pairs=[];
   for(let j=1;j<rs.length;j++){const a=rs[j-1],b=rs[j];if(parseInt(b.r,10)-parseInt(a.r,10)===1&&points[j-1].distanceTo(points[j])<5)pairs.push([points[j-1],points[j]]);}
   const color=palette[i%palette.length],trace=new T.Group();trace.userData.chain=chain;this.backbone.add(trace);let run=[];const flush=()=>{if(run.length>1){const mesh=new T.Mesh(new T.TubeGeometry(new T.CatmullRomCurve3(run,false,'centripetal'),Math.max(8,run.length*3),.45,8,false),material(color));mesh.userData={flowProtein:true,chain,trace:true};trace.add(mesh);}run=[];};
   for(let j=0;j<points.length;j++){if(j&&(parseInt(rs[j].r,10)-parseInt(rs[j-1].r,10)!==1||points[j].distanceTo(points[j-1])>=5))flush();run.push(points[j]);}flush();
   const dots=instances(new T.SphereGeometry(1,10,8),color,points,this.backbone,.3);dots.userData={flowProtein:true,chain,offset};offset+=rs.length;this.pickMeshes.push(dots);
  }
  this.atoms=this.pickMeshes[0];this.sourcePoints=rows.map(r=>r.p);this.updateDisplayPoints();
  const ring=new T.LineLoop(new T.BufferGeometry().setFromPoints(Array.from({length:65},(_,i)=>new T.Vector3(Math.cos(i/64*Math.PI*2)*this.sourceExtent*.65,Math.sin(i/64*Math.PI*2)*this.sourceExtent*.65,-this.sourceExtent*.54))),new T.LineBasicMaterial({color:'#a88659',transparent:true,opacity:.5}));ring.quaternion.copy(this.protein.quaternion).invert();this.protein.add(ring);
  const ruler=new T.Line(new T.BufferGeometry().setFromPoints([new T.Vector3(-10,-this.sourceExtent*.65,0),new T.Vector3(10,-this.sourceExtent*.65,0)]),new T.LineBasicMaterial({color:'#fff3d5'}));ruler.quaternion.copy(this.protein.quaternion).invert();this.protein.add(ruler);
  this.highlightMesh=new T.Mesh(new T.SphereGeometry(.85,12,8),material('#ffffff'));this.highlightMesh.visible=false;this.protein.add(this.highlightMesh);
 }
 updateDisplayPoints(){this.displayPoints=this.sourcePoints.map(p=>new T.Vector3(...p).sub(this.center).applyQuaternion(this.protein.quaternion).multiplyScalar(this.factor).add(this.protein.position).toArray());}
 setRepresentation(value){
  if(!this.record)return;this.representation=value;if(value==='atoms'){this.selectedChain='all';this.backbone.children.forEach(o=>o.visible=true);}
  if(value==='atoms'&&!this.atomic){
   this.atomic=new T.Group();this.protein.add(this.atomic);const points=this.rawAtoms.map(a=>new T.Vector3(a.x,a.y,a.z).sub(this.center));
   const dots=instances(new T.SphereGeometry(1,8,6),'#b5d5d1',points,this.atomic,.4),color=new T.Color();const colors={C:'#b7c7bd',N:'#8ab4ec',O:'#ec9d91',S:'#edcf7f',P:'#d7aaee'};
   this.rawAtoms.forEach((a,i)=>dots.setColorAt(i,color.set(colors[a.elem]||'#bbb7b1')));dots.userData={flowProtein:true,atomic:true};this.atomicMesh=dots;
   const pairs=[];this.rawAtoms.forEach((a,i)=>{for(const j of a.bonds||[])if(j>i&&points[j]&&points[i].distanceTo(points[j])<2.5)pairs.push([points[i],points[j]]);});cylinders(pairs,.11,'#9dafb0',this.atomic);
  }
  this.backbone.visible=value!=='atoms';if(this.atomic)this.atomic.visible=value==='atoms';
 }
 setChain(chain){if(!this.record)return;if(chain!=='all'&&this.representation==='atoms')this.setRepresentation('backbone');this.selectedChain=chain;this.backbone.children.forEach((o,i)=>o.visible=chain==='all'||Object.keys(this.record.chains)[Math.floor(i/2)]===chain);if(this.atomic)this.atomic.visible=this.representation==='atoms'&&chain==='all';if(chain!=='all'&&this.representation==='atoms'){this.backbone.visible=true;}}
 selectResidue(chain,resid){const r=this.residues?.find(r=>r.chain===chain&&String(r.r)===String(resid));if(!r)return;this.highlightMesh.position.copy(new T.Vector3(...r.p).sub(this.center));this.highlightMesh.visible=true;return r;}
 pick(ray){const visible=this.representation==='atoms'&&this.atomic?.visible?[this.atomicMesh]:this.pickMeshes.filter(m=>m.visible);const hit=ray.intersectObjects(visible,false)[0];if(!hit){const traceHit=ray.intersectObject(this.backbone,true).find(h=>h.object.userData.trace&&(!this.selectedChain||this.selectedChain==='all'||h.object.userData.chain===this.selectedChain));if(!traceHit)return null;const p=this.protein.worldToLocal(traceHit.point.clone()).add(this.center);let nearest=null,d=Infinity;for(const residue of this.residues){if(residue.chain!==traceHit.object.userData.chain)continue;const x=p.distanceToSquared(new T.Vector3(...residue.p));if(x<d){d=x;nearest=residue;}}if(nearest)this.selectResidue(nearest.chain,nearest.r);return nearest;}let r;if(hit.object.userData.atomic){const a=this.rawAtoms[hit.instanceId];r=this.residues.find(r=>r.chain===a.chain&&parseInt(r.r,10)===a.resi);}else r=this.residues[hit.object.userData.offset+hit.instanceId];if(r)this.selectResidue(r.chain,r.r);return r||{chain:'?',r:'?',n:'atom'};}
 setLevel(value){
  this.active=true;this.targetLevel=T.MathUtils.clamp(+value,0,2);
  const x=this.targetLevel,a=T.MathUtils.smoothstep(x,0,1),b=T.MathUtils.smoothstep(x,1,2);
  const p=this.root.localToWorld(this.protein.position.clone()),patch=this.anchor.clone().lerp(p,.25),target=this.startTarget.clone().lerp(patch,a).lerp(p,b);
  const side=new T.Vector3(1,0,0).applyQuaternion(this.root.getWorldQuaternion(new T.Quaternion()));
  const patchDistance=this.patchRadius*(this.camera.aspect<1?7.8:5.4),proteinDistance=this.patchRadius*(this.camera.aspect<1?3.6:2.65);
  const near=patch.clone().addScaledVector(this.normal,patchDistance).addScaledVector(side,patchDistance*.12),close=p.clone().addScaledVector(this.normal,proteinDistance).addScaledVector(side,proteinDistance*.12);
  const position=this.start.clone().lerp(near,a).lerp(close,b);
  this.transition={from:this.camera.position.clone(),targetFrom:this.controls.target.clone(),position,target,levelFrom:this.level,elapsed:0,last:performance.now()};
 }
 stop(){this.transition=null;this.active=true;this.controls.enabled=true;}
 update(){
  const transition=this.transition;if(transition){const now=performance.now();transition.elapsed+=(now-transition.last)/1000;transition.last=now;const duration=(matchMedia('(prefers-reduced-motion: reduce)').matches||document.documentElement.dataset.reduceMotion==='true')?0:1.05,t=duration?Math.min(1,transition.elapsed/duration):1,k=t*t*(3-2*t);
   this.camera.position.lerpVectors(transition.from,transition.position,k);this.controls.target.lerpVectors(transition.targetFrom,transition.target,k);this.level=T.MathUtils.lerp(transition.levelFrom,this.targetLevel,k);if(t===1)this.transition=null;
  }
  this.camera.near=Math.max(.00001,this.camera.position.distanceTo(this.controls.target)*.001);this.camera.far=Math.max(100,this.start.length()*5);this.camera.updateProjectionMatrix();this.controls.update();
  this.contextMeshes.forEach(o=>o.visible=this.level<.5);
  this.bodySurfaces.forEach(o=>{o.material.transparent=true;o.material.opacity=this.level>.5?.88:.55;o.material.depthWrite=this.level>.5;});
  this.onStatus?.(this.level);
 }
 info(){
  let maxPairError=0;if(this.sourcePoints.length>1){for(let i=0;i<Math.min(64,this.sourcePoints.length-1);i++){const j=this.sourcePoints.length-1-i;const expected=new T.Vector3(...this.sourcePoints[i]).distanceTo(new T.Vector3(...this.sourcePoints[j]))*this.factor,actual=new T.Vector3(...this.displayPoints[i]).distanceTo(new T.Vector3(...this.displayPoints[j]));maxPairError=Math.max(maxPairError,Math.abs(expected-actual));}}
  let renderedCoordinateError=0;const m=new T.Matrix4();for(const mesh of this.pickMeshes||[]){for(const i of [0,Math.floor(mesh.count/2),mesh.count-1]){if(i<0)continue;mesh.getMatrixAt(i,m);const expected=new T.Vector3(...this.sourcePoints[mesh.userData.offset+i]).sub(this.center),actual=new T.Vector3().setFromMatrixPosition(m);renderedCoordinateError=Math.max(renderedCoordinateError,expected.distanceTo(actual));}}
  let surfaceResidual=0;for(const sample of this.patchSamples){const expected=new T.Vector3();sample.triangle.forEach((id,k)=>expected.addScaledVector(this.graph.vertices[id].p,sample.weights[k]));const actual=this.root.localToWorld(sample.p.clone());surfaceResidual=Math.max(surfaceResidual,expected.distanceTo(actual));}
  return {surfaceResidual,renderedCoordinateError,scene:this.group.uuid,region:this.root.uuid,proteinObject:this.protein.uuid,regionSerial:this.regionSerial,active:this.active,level:this.level,target:this.targetLevel,transitioning:!!this.transition,source:this.record?.id||null,sha256:this.record?.sha256||null,chains:Object.keys(this.record?.chains||{}),residues:this.sourcePoints.length,atoms:this.rawAtoms.length,uniformScale:this.factor||null,anchor:this.anchor.toArray(),localRadius:this.localRadius,patchRadius:this.patchRadius,surfaceSamples:this.patchSamples.length,surfaceGeometry:this.surface.geometry.type,boundary:this.surface.userData.boundaryLabel,representation:this.representation,placement:'Linked inspection cradle; NO native insertion, localization or pose',coordinatePreservationMaxPairError:maxPairError,camera:this.camera.position.toArray(),targetPoint:this.controls.target.toArray(),near:this.camera.near};
 }
}
