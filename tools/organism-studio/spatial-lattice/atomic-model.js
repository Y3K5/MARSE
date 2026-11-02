// Read-only optional atom context. Original producer pins and residue readers remain authoritative.
import {residueRows} from './residue-model.js';
export const atomicKinds={component_template:'Component topology',polymer_link:'Polymer sequence + geometry',component_coordination:'Component coordination',source_covale:'Source-declared covalent link',source_disulf:'Source-declared disulfide',source_metalc:'Source-declared coordination',source_hydrog:'Source-declared hydrogen connection'};
export const atomicOrders={sing:'Single',doub:'Double',trip:'Triple',quad:'Quadruple',arom:'Aromatic',delo:'Delocalized',pi:'Pi',poly:'Polymeric',unknown:'Order unassigned'};
export const elementColors={C:'#aab8bf',N:'#5799ef',O:'#f06e6b',S:'#e9ce62',P:'#e59b55',H:'#ebefed',FE:'#d98b5c',ZN:'#91baca',CA:'#9bbf91',MG:'#80cda4'};
const sha=async raw=>Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',raw))).map(n=>n.toString(16).padStart(2,'0')).join('');
export function parseAtomicPdb(text){
 const atoms=[],seen=new Set();for(const line of text.split(/\r?\n/)){
  if(line.startsWith('ENDMDL'))break;if(!/^(ATOM  |HETATM)/.test(line))continue;
  const altloc=line.slice(16,17).trim(),occupancy=line.slice(54,60).trim()?Number(line.slice(54,60)):null;if(!['','A'].includes(altloc)||occupancy!==null&&occupancy<=0)continue;
  const atom={i:atoms.length,serial:Number(line.slice(6,11)),name:line.slice(12,16).trim(),element:line.slice(76,78).trim().toUpperCase(),chain:line.slice(21,22).trim()||'_',residue:Number(line.slice(22,26)),insertion:line.slice(26,27).trim(),resname:line.slice(17,20).trim(),altloc,xyz:[Number(line.slice(30,38)),Number(line.slice(38,46)),Number(line.slice(46,54))],occupancy,b:line.slice(60,66).trim()?Number(line.slice(60,66)):null,group:line.slice(0,6).trim()};
  const required=[line.slice(6,11),line.slice(22,26),line.slice(30,38),line.slice(38,46),line.slice(46,54)];const key=JSON.stringify([atom.chain,atom.residue,atom.insertion,atom.resname,atom.name]);if(required.some(v=>!v.trim())||!atom.name||!atom.resname||[atom.occupancy,atom.b].some(v=>v!==null&&!Number.isFinite(v))||seen.has(key)||!Number.isInteger(atom.serial)||!Number.isInteger(atom.residue)||!atom.element||atom.xyz.some(x=>!Number.isFinite(x)))throw Error('Atomic source identities or coordinates are incomplete.');seen.add(key);atoms.push(atom);
 }if(!atoms.length)throw Error('No primary retained atoms.');return atoms;
}
export function atomicResidueKey(a){return JSON.stringify([a.chain,a.residue,a.insertion,a.resname]);}
export function atomicLabel(a){return `${a.name} (${a.element}) · ${a.resname} ${a.chain}:${a.residue}${a.insertion}`;}
export function validateAtomicPayload(p,pdbText,data,ledger,id){
 const shape=data.shape_traces[id];if(p.schema!=='marse-atomic-context-v1'||p.protein_id!==id||!shape||p.source_view_sha256!==shape.source_view_sha256||p.scientific_execution!==false||!Array.isArray(p.atoms)||!Array.isArray(p.bonds))throw Error('Atomic model differs from this exact retained fold.');
 const url=new URL(p.context?.source_url);if(!['http:','https:'].includes(url.protocol)||url.username||url.password||!['experimental','prediction'].includes(p.context.kind))throw Error('Atomic source context is unsafe or unspecified.');
 const source=parseAtomicPdb(pdbText),keys=['i','serial','name','element','chain','residue','insertion','resname','altloc','xyz','occupancy','b','group'];if(source.length!==p.atoms.length)throw Error('Atomic source coverage differs.');
 for(const [i,a]of p.atoms.entries())for(const k of keys)if(JSON.stringify(a[k])!==JSON.stringify(source[i][k]))throw Error('Atomic source name/element/coordinates differ.');
 const anchors=source.filter(a=>a.group==='ATOM'&&a.name==='CA'),rows=residueRows(shape,ledger.entries[id]);if(anchors.length!==rows.length)throw Error('Atomic/backbone coverage differs.');
 const by=new Map(anchors.map(a=>[JSON.stringify([a.chain,a.residue,a.insertion]),a]));for(const r of rows){const a=by.get(JSON.stringify([r.chain,r.number,r.insertion]));if(!a||a.resname!==r.name||Math.hypot(...a.xyz.map((v,k)=>v-r.point[k+1]))>.00001)throw Error('Atomic/backbone correspondence held.');}
 const seen=new Set();for(const b of p.bonds){if(!Number.isInteger(b.a)||!Number.isInteger(b.b)||b.a<0||b.b>=p.atoms.length||b.a>=b.b||!Object.hasOwn(atomicKinds,b.kind)||!['sing','doub','trip','quad','arom','delo','pi','poly','unknown'].includes(b.order))throw Error('Unsupported atomic connection.');const key=b.a+':'+b.b;if(seen.has(key))throw Error('Duplicate atomic connection.');seen.add(key);const d=Math.hypot(...p.atoms[b.a].xyz.map((v,k)=>v-p.atoms[b.b].xyz[k]));if(!Number.isFinite(b.distance)||Math.abs(d-b.distance)>.0001)throw Error('Atomic connection geometry differs.');}
 return p;
}
export function atomicNeighbors(payload,index){return payload.bonds.filter(b=>b.a===index||b.b===index).map(b=>({atom:payload.atoms[b.a===index?b.b:b.a],bond:b}));}
export function atomicEnvironment(payload,index,radius=6,chain='all'){
 const selected=payload.atoms[index];if(!selected||!Number.isFinite(radius)||radius<0||radius>20)throw Error('Invalid atom environment.');return payload.atoms.filter(a=>(chain==='all'||a.chain===chain)&&Math.hypot(...a.xyz.map((v,k)=>v-selected.xyz[k]))<=radius).map(a=>a.i);
}
export class AtomicStore{
 constructor(data,ledger,sourceHash,read=url=>fetch(url),residueHash=null){Object.assign(this,{data,ledger,sourceHash,residueHash,read,cache:new Map()});}
 async bytes(path){const r=await this.read('runtime/atomic/'+path);if(!r.ok)throw Error('No optional local atomic import is available.');const bytes=await r.arrayBuffer();if(bytes.byteLength>64*1024*1024)throw Error('Atomic payload exceeds size limit.');return bytes;}
 async index(){if(this.catalog)return this.catalog;const receipt=JSON.parse(new TextDecoder().decode(await this.bytes('IMPORT.json'))),raw=await this.bytes('index.json');if(receipt.schema!=='marse-atomic-import-v1'||await sha(raw)!==receipt.index_sha256)throw Error('Atomic import index checksum differs.');const index=JSON.parse(new TextDecoder().decode(raw));if(index.schema!=='marse-atomic-index-v1'||index.source_snapshot_sha256!==this.sourceHash||this.residueHash&&index.residue_ledger_sha256!==this.residueHash||index.scientific_execution!==false)throw Error('Atomic import belongs to another source snapshot.');this.receipt=receipt;this.catalog=index;return index;}
 async load(id){
  if(this.cache.has(id))return this.cache.get(id);const index=await this.index(),entry=index.entries[id],shape=this.data.shape_traces[id];if(!entry||!shape||entry.source_view_sha256!==shape.source_view_sha256)throw Error('No full-atom import for this exact saved fold. No related protein is substituted.');
  if(!/^[a-f0-9]{64}\.json$/.test(entry.file)||entry.pdb_file!==shape.source_view_sha256+'.pdb')throw Error('Unsafe atomic payload identity.');const [raw,pdb]=await Promise.all([this.bytes(entry.file),this.bytes(entry.pdb_file)]);if(await sha(raw)!==entry.sha256||this.receipt.files[entry.file]!==entry.sha256||await sha(pdb)!==shape.source_view_sha256||this.receipt.files[entry.pdb_file]!==shape.source_view_sha256)throw Error('Atomic source or payload checksum differs.');
  const p=validateAtomicPayload(JSON.parse(new TextDecoder().decode(raw)),new TextDecoder().decode(pdb),this.data,this.ledger,id);this.cache.set(id,p);return p;
 }
}
