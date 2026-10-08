import {readFile,writeFile,mkdir} from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';
import {fileURLToPath} from 'node:url';
const root=fileURLToPath(new URL('.',import.meta.url));process.chdir(root);
// Read-only teaching endpoints; the full historical trajectory is held outside Git.
const input=JSON.parse(await readFile('src/pocket-input.json','utf8'));
// Deduplicate shared anatomy and pack viewer buffers; canonical Blender JSON stays intact.
const authored=JSON.parse(await readFile('assets/pocket-organic/anatomy.json','utf8'));
const pool=[],lookup=new Map(),cases={};
for(const [health,parts] of Object.entries(authored.cases))cases[health]=parts.map(p=>{
 const key=createHash('sha256').update(JSON.stringify(p)).digest('hex');if(lookup.has(key))return lookup.get(key);
 const index=pool.length;lookup.set(key,index);pool.push({id:p.id,name:p.name,positions:Buffer.from(new Float32Array(p.positions).buffer).toString('base64'),indices:Buffer.from(new Uint16Array(p.indices).buffer).toString('base64')});
 if(p.indices.some(i=>!Number.isInteger(i)||i<0||i>65535))throw Error('Mesh exceeds Uint16 budget');return index;
});
await writeFile('src/pocket-organic-input.json',JSON.stringify({format:'marse.pocket-organic-packed/1',units:authored.units,calibration:null,precision:'float32 viewer vertices, Uint16 indices; little-endian',source_sha256:createHash('sha256').update(await readFile('assets/pocket-organic/anatomy.json')).digest('hex'),pool,cases})+'\n');
await mkdir('dist',{recursive:true});
await build({entryPoints:['src/core.ts'],outfile:'dist/core.mjs',bundle:true,platform:'node',format:'esm',target:'es2022',external:['three']});
await build({entryPoints:['src/pocket.ts'],outfile:'dist/pocket-explorer.js',bundle:true,minify:true,format:'iife',target:'es2022',legalComments:'eof'});
await build({entryPoints:['src/pocket-geometry.ts','src/pocket-science.ts','src/pocket-section.ts'],outdir:'dist/pocket-core',outExtension:{'.js':'.mjs'},bundle:true,platform:'node',format:'esm',target:'es2022',external:['three']});
const files=['pocket.html','pocket-explorer.css','dist/pocket-explorer.js','src/pocket.ts','src/pocket-geometry.ts','src/pocket-science.ts','src/pocket-input.json','src/pocket-materials.ts','src/pocket-section.ts','src/pocket-organic-input.json','assets/pocket-organic/anatomy.json','tools/build-organic-pocket.py','build.mjs','package.json','package-lock.json'];
const assets=[];for(const file of files){const b=await readFile(file);assets.push({file,bytes:b.length,sha256:createHash('sha256').update(b).digest('hex')});}
await writeFile('dist/pocket-explorer-manifest.json',JSON.stringify({format:'marse.pocket-explorer-manifest/1',version:'0.1.0',host:'human',tooth:'FDI 36',geometry:'Authored, uncalibrated; expert review pending',units:'illustrative model units',source:input.source,files:assets},null,2)+'\n');
console.log('Pocket Explorer built locally. Two exact saved endpoints, nine taxa. Studio bundle untouched.');
