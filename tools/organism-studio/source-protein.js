// Reuse only the locally vendored molecular parser; it never creates a viewer.
const cache=new Map();
export async function loadSourceProtein(id){
 if(cache.has(id))return cache.get(id);
 const response=await fetch('structures/'+id+'.json');if(!response.ok)throw Error(id+' metadata unavailable');
 const record=await response.json(),rawResponse=await fetch('structures/'+id+'.'+(record.format||'pdb'));
 if(!rawResponse.ok)throw Error(id+' coordinate file unavailable');const bytes=await rawResponse.arrayBuffer();
 const digest=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes))).map(x=>x.toString(16).padStart(2,'0')).join('');
 if(digest!==record.sha256)throw Error(id+' coordinate checksum mismatch; source not displayed');
 const parser=new $3Dmol.GLModel(0,{});parser.addMolData(new TextDecoder().decode(bytes),record.format||'pdb',{keepH:true,doAssembly:false});
 const allowed=new Set(Object.keys(record.chains)),all=parser.selectedAtoms({}),retained=all.filter(a=>allowed.has(a.chain)&&!['HOH','WAT'].includes(a.resn)),remap=new Map(retained.map((a,i)=>[a.index,i]));
 const atoms=retained.map(a=>({...a,bonds:(a.bonds||[]).filter(i=>remap.has(i)).map(i=>remap.get(i))}));
 if(!atoms.length)throw Error(id+' has no atoms in the recorded source chains');
 const result={record,atoms};cache.set(id,result);return result;
}
