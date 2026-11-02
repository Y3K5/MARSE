/* Optional private extensions for the organism studio.
   Reads private/extensions.json, which exists only when serve.py runs with --private DIR (a folder kept outside
   this repository). Without it this script only restores a deep-linked selection and does nothing else.
   extensions.json: {"links": [{"label": "...", "href": "private/..."}], "structuralPieces": "private/.../studio.html"}
   Only same-origin hrefs under private/ are accepted. */
(()=>{
  const safe=h=>typeof h==='string'&&/^private\/[A-Za-z0-9_./%-]+$/.test(h)&&!h.split('/').includes('..');

  // Deep links (?organism=&component=) restore a selection once the organism list has loaded.
  const query=new URLSearchParams(location.search),organism=query.get('organism'),component=query.get('component');
  if(organism){
    let applied=false;const select=document.getElementById('organism');
    const restore=()=>{if(applied||![...select.options].some(o=>o.value===organism))return;applied=true;select.value=organism;select.dispatchEvent(new Event('change',{bubbles:true}));const protein=document.getElementById('protein');if(component!==null&&[...protein.options].some(o=>o.value===component)){protein.value=component;protein.dispatchEvent(new Event('change',{bubbles:true}));}};
    new MutationObserver(restore).observe(select,{childList:true});restore();
  }

  // Structural hand-off: opens the selected source structure in a private Structural Pieces workspace.
  function mountHandoff(studio){
    const map=document.querySelector('.connection-map');if(!map)return;
    const section=document.createElement('section');section.className='structural-handoff';
    const link=document.createElement('a');link.id='open-structural-piece';link.href=studio;link.textContent='View selected structural piece';
    const note=document.createElement('p');note.id='structural-handoff-note';note.textContent='Retained molecular sources open in the private Structural Pieces workspace.';
    section.append(link,note);map.append(section);
    function sync(){
      const org=document.getElementById('organism').value;let source=null;
      try{const url=new URL(document.querySelector('#structure-info a')?.href||'',location.href);if(/(^|\.)rcsb\.org$/.test(url.hostname))source=url.pathname.match(/^\/structure\/([0-9a-z]{4})\/?$/i)?.[1].toUpperCase();}catch{}
      if(!org||!source){link.href=studio;link.textContent='Browse structural pieces';note.textContent='No source coordinates loaded for this selection. Browse the native catalogue; missing components remain explicit.';return;}
      link.href=studio+'?'+new URLSearchParams({source,organism:org,component:document.getElementById('protein').value});
      link.textContent='View '+source+' in Structural Pieces';note.textContent='Same retained coordinate source; molecular size is shown in Å. Native pose remains unknown.';
    }
    link.addEventListener('click',sync);
    const info=document.getElementById('structure-info');if(info)new MutationObserver(sync).observe(info,{childList:true,subtree:true});
    for(const id of ['protein','organism'])document.getElementById(id)?.addEventListener('change',()=>queueMicrotask(sync));
    sync();
  }

  fetch('private/extensions.json',{cache:'no-store'}).then(r=>r.ok?r.json():null).then(ext=>{
    if(!ext||typeof ext!=='object')return;
    const nav=document.createElement('nav');nav.className='workspace-nav';nav.setAttribute('aria-label','Private workspaces');
    for(const l of Array.isArray(ext.links)?ext.links:[]){if(!safe(l&&l.href)||typeof l.label!=='string')continue;const a=document.createElement('a');a.href=l.href;a.textContent=l.label;nav.append(a);}
    if(nav.children.length)document.querySelector('header .local')?.before(nav);
    if(safe(ext.structuralPieces))mountHandoff(ext.structuralPieces);
  }).catch(()=>{});
})();
