// Local, texture-free GLB export of the same qualitative bodies used in the viewer.
import {mkdir,writeFile} from 'node:fs/promises';
import {resolve,join} from 'node:path';
import {createHash} from 'node:crypto';
import * as T from 'three';
import {microbialObject} from '../dist/pocket-core/pocket-crevice.mjs';
import {MORPHOLOGY} from '../dist/pocket-core/pocket-morphology.mjs';
import {TAXA} from '../dist/pocket-core/pocket-science.mjs';
const output=resolve(process.argv[2]||'assets/crevice-objects');await mkdir(output,{recursive:true});
function glb(t){
 const g=microbialObject(t.id),chunks=[],views=[],accessors=[];let offset=0;
 const add=(array,componentType,type,target)=>{const raw=Buffer.from(array.buffer,array.byteOffset,array.byteLength),pad=Buffer.alloc((4-raw.length%4)%4);const view=views.length;views.push({buffer:0,byteOffset:offset,byteLength:raw.length,target});chunks.push(raw,pad);offset+=raw.length+pad.length;const a={bufferView:view,componentType,count:array.length/(type==='VEC3'?3:1),type};accessors.push(a);return accessors.length-1;};
 const p=g.getAttribute('position'),n=g.getAttribute('normal'),pos=add(new Float32Array(p.array),5126,'VEC3',34962),normal=add(new Float32Array(n.array),5126,'VEC3',34962);
 g.computeBoundingBox();accessors[pos].min=g.boundingBox.min.toArray();accessors[pos].max=g.boundingBox.max.toArray();
 const primitive={attributes:{POSITION:pos,NORMAL:normal},material:0,mode:4};if(g.index)primitive.indices=add(new Uint32Array(g.index.array),5125,'SCALAR',34963);
 const c=new T.Color(t.color),m=MORPHOLOGY[t.id];
 const document={asset:{version:'2.0',generator:'MARSE qualitative microbial-object exporter'},scene:0,scenes:[{nodes:[0]}],nodes:[{name:t.name,mesh:0,extras:{taxon:t.id,identificationLevel:t.level,source:m.source,context:m.context,claim:'Qualitative illustrative morphology; artistic color and relief; no physical scale, contacts, molecular sites or inferred location'}}],meshes:[{name:t.name,primitives:[primitive]}],materials:[{name:t.short,pbrMetallicRoughness:{baseColorFactor:[c.r,c.g,c.b,1],metallicFactor:0,roughnessFactor:.38},doubleSided:true}],buffers:[{byteLength:offset}],bufferViews:views,accessors,extras:{units:'normalized authoring units, not physical dimensions',coordinateSystem:'right-handed Y-up; body centered at origin',savedFieldsUsed:false}};
 let json=Buffer.from(JSON.stringify(document));json=Buffer.concat([json,Buffer.alloc((4-json.length%4)%4,32)]);const binary=Buffer.concat(chunks),header=Buffer.alloc(12),jHeader=Buffer.alloc(8),bHeader=Buffer.alloc(8);header.writeUInt32LE(0x46546c67);header.writeUInt32LE(2,4);header.writeUInt32LE(12+16+json.length+binary.length,8);jHeader.writeUInt32LE(json.length);jHeader.writeUInt32LE(0x4e4f534a,4);bHeader.writeUInt32LE(binary.length);bHeader.writeUInt32LE(0x004e4942,4);g.dispose();return Buffer.concat([header,jHeader,json,bHeader,binary]);
}
const files=[];for(const t of TAXA){const bytes=glb(t),name=t.id+'.glb';await writeFile(join(output,name),bytes);files.push({file:name,taxon:t.id,name:t.name,source:MORPHOLOGY[t.id].source,bytes:bytes.length,sha256:createHash('sha256').update(bytes).digest('hex')});}
await writeFile(join(output,'objects-manifest.json'),JSON.stringify({format:'marse.illustrative-microbial-objects/1',claim:'Qualitative shape examples; not measured cells or colonization.',physicalScale:null,savedFieldsUsed:false,files},null,2)+'\n');console.log('Exported nine local illustrative GLB objects.');
