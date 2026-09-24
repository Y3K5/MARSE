"""Self-contained HTML viewer for exported ecosystem frames.

The embedded JavaScript is intentionally kept in one standalone document so
the result can be opened without a web server or frontend dependency.
"""

# The generated JavaScript is formatted for the browser, not as Python.
# ruff: noqa: E501

from __future__ import annotations

from pathlib import Path

from marse.ecosystem.model import EcosystemResult


def write_viewer(result: EcosystemResult, path: str | Path) -> Path:
    """Write a dependency-free browser viewer with embedded simulation data."""
    viewer_path = Path(path)
    data = result.write_frames(viewer_path.with_name("frames.json"))
    html = f"""<!doctype html>
<meta charset="utf-8"><title>MARSE ecosystem viewer</title>
<style>
:root{{color-scheme:dark;--panel:#171d27;--line:#334155;--text:#e7edf5;--muted:#9fb0c4}}
*{{box-sizing:border-box}} body{{font:15px system-ui;background:radial-gradient(circle at top,#172334,#0b0f15 70%);color:var(--text);margin:0;padding:24px}}
main{{max-width:1200px;margin:auto}} h1{{font-size:clamp(1.3rem,3vw,2rem);margin:0 0 4px}}
.subtitle,.mut{{color:var(--muted)}} .panel{{background:color-mix(in srgb,var(--panel) 90%,transparent);border:1px solid var(--line);border-radius:14px;padding:14px;margin-top:14px;box-shadow:0 12px 35px #0004}}
#controls{{display:flex;flex-wrap:wrap;align-items:center;gap:8px}} button,input,select{{accent-color:#46b3ff;background:#202a38;color:var(--text);border:1px solid #40516a;border-radius:7px;padding:6px 9px}}
button{{cursor:pointer}} button:hover{{border-color:#7dd3fc}} label{{color:var(--muted)}} #canvas-wrap{{overflow:auto;max-width:100%;border-radius:10px;background:#080b10}} canvas{{border:1px solid #526071;border-radius:10px;cursor:crosshair;background:#080b10;display:block;image-rendering:auto}}
.layout{{display:grid;grid-template-columns:minmax(0,1fr) 280px;gap:14px;align-items:start}} #legend{{font-family:ui-monospace,monospace;color:#b8c7d9}} #inspect{{min-height:3em;color:#b8c7d9;line-height:1.5}}
#stats{{line-height:1.6;color:#d8e3ef}} .metric{{display:inline-block;margin-right:14px}} .swatch{{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:5px}}
@media(max-width:800px){{body{{padding:12px}}.layout{{grid-template-columns:1fr}}canvas{{width:100%}}}}
</style>
<main><h1>MARSE ecosystem: {result.config.experiment_id}</h1>
<div class="subtitle">Interactive spatial simulation explorer · <span id="status">paused</span></div>
<div class="panel" id="controls"><button id="play">Play</button><button id="back">Back</button><button id="forward">Next</button><label>Speed <select id="speed"><option value="240">0.25x</option><option value="100" selected>1x</option><option value="35">3x</option><option value="10">10x</option></select></label><label>View <select id="render-mode"><option value="combined" selected>Fields + particles</option><option value="smooth">Smooth fields</option><option value="cells">Cell grid</option><option value="agents">Agents only</option></select></label><label>Zoom <input id="zoom" type="range" min="1" max="8" step="0.5" value="3"></label><label>Opacity <input id="opacity" type="range" min="0.25" max="1" step="0.05" value="1"></label><input id="time" type="range" min="0" max="{len(result.frames) - 1}" value="0" style="flex:1;min-width:180px">
<span id="label"></span><label>Layer <select id="layer"></select></label><span id="legend"></span></div>
<div class="layout"><section class="panel"><div id="stats"></div><div id="canvas-wrap"><canvas id="view" width="{result.config.width}" height="{result.config.height}"></canvas></div></section>
<aside class="panel"><strong>Cell inspector</strong><div id="inspect">Click a cell to inspect its local state.</div><hr><strong>Layers</strong><p class="mut">Species, nutrients, conditions, additives, quorum phenotypes, effective growth rates, and limiting factors.</p><p class="mut">Keyboard: Space play/pause · ←/→ step · Home/End jump.</p></aside></div>
<p class="mut">Particles are deterministic population representatives sampled from biomass density; explicit immune agents are shown when supplied by a frame. Yellow outlines mark mutations; green outlines mark positive rates.</p></main>
<script>
const data = {data.read_text(encoding="utf-8")};
const canvas=document.getElementById('view'), ctx=canvas.getContext('2d'), stats=document.getElementById('stats'), inspect=document.getElementById('inspect'), status=document.getElementById('status');
const slider=document.getElementById('time'), label=document.getElementById('label'), selector=document.getElementById('layer'), legend=document.getElementById('legend'), speed=document.getElementById('speed'), renderMode=document.getElementById('render-mode'), zoom=document.getElementById('zoom'), opacity=document.getElementById('opacity');
let playing=false, timer;
const colors=['#46b3ff','#ff6b6b','#7ee081','#c084fc','#fb923c','#f472b6'];
const names=[...data.species], nutrientNames=[...data.nutrients], conditionNames=[...(data.conditions||[])], additiveNames=[...(data.additives||[])];
const layers=[];
for(const [kind, values] of [['species',names],['nutrient',nutrientNames],['condition',conditionNames],['additive',additiveNames]]) for(const [i,name] of values.entries()) layers.push({{kind,index:i,label:`${{name}} (${{kind}})`}});
for(const [i,name] of names.entries()) layers.push({{kind:'phenotype',index:i,label:`${{name}} (phenotype state)`}});
for(const [i,name] of names.entries()) layers.push({{kind:'rate',index:i,label:`${{name}} (effective growth rate)`}});
for(const [i,name] of names.entries()) layers.push({{kind:'limit',index:i,label:`${{name}} (limiting factor)`}});
for(const [i,layer] of layers.entries()) selector.innerHTML+=`<option value="${{i}}">${{layer.label}}</option>`;
function selectedLayer(){{return layers[+selector.value]}}
function fieldFor(frame, layer){{
 if(layer.kind==='species') return frame.species[names[layer.index]];
 if(layer.kind==='nutrient') return frame.nutrients[nutrientNames[layer.index]];
 if(layer.kind==='condition') return frame.conditions[conditionNames[layer.index]];
 if(layer.kind==='additive') return frame.additives[additiveNames[layer.index]];
 if(layer.kind==='phenotype') return frame.phenotypes[layer.index].map(row=>row.map(value=>`state ${{value}}`));
 if(layer.kind==='rate') return frame.niche.effective_growth_rate_per_h[names[layer.index]];
 return frame.niche.limiting_factor[names[layer.index]];
}}
function color(value, layer, max){{
 if(layer.kind==='limit'){{const palette=['#ef4444','#f59e0b','#22c55e','#38bdf8','#c084fc']; return palette[Math.abs(value.split('').reduce((a,c)=>a+c.charCodeAt(0),0))%palette.length]}}
 if(layer.kind==='phenotype'){{const state=Number(String(value).replace('state ','')); return ['#334155','#f59e0b','#ec4899','#22c55e','#a78bfa'][state%5]}}
 const ratio=Math.max(0,Math.min(1,value/Math.max(max,1e-12))); return `rgb(${{Math.round(40+180*ratio)}},${{Math.round(90+130*ratio)}},${{Math.round(150+90*ratio)}})`;
}}
function drawField(field,layer,max){{
 const smooth=renderMode.value==='smooth'||renderMode.value==='combined', width=data.width, height=data.height;
 ctx.clearRect(0,0,width,height); ctx.imageSmoothingEnabled=smooth; ctx.globalAlpha=+opacity.value;
 if(smooth){{
  const source=document.createElement('canvas'); source.width=width; source.height=height;
  const sourceCtx=source.getContext('2d'); sourceCtx.imageSmoothingEnabled=false;
  for(let y=0;y<height;y++)for(let x=0;x<width;x++){{sourceCtx.fillStyle=color(field[y][x],layer,max);sourceCtx.fillRect(x,y,1,1);}}
  ctx.drawImage(source,0,0,width,height);
 }} else {{
  for(let y=0;y<height;y++)for(let x=0;x<width;x++){{ctx.fillStyle=color(field[y][x],layer,max);ctx.fillRect(x,y,1,1);}}
 }} ctx.globalAlpha=1;
}}
function seeded(seed){{const value=Math.sin(seed*12.9898)*43758.5453;return value-Math.floor(value)}}
function populationParticles(frame){{
 const particles=[]; const maxParticles=700;
 for(let species=0;species<names.length;species++){{
  const field=frame.species[names[species]], occupied=[];
  for(let y=0;y<data.height;y++)for(let x=0;x<data.width;x++)if(field[y][x]>0.01)occupied.push([y,x,field[y][x]]);
  const stride=Math.max(1,Math.ceil(occupied.length*4/maxParticles));
  occupied.filter((_,index)=>index%stride===0).forEach(([y,x,density],index)=>{{
   const count=Math.max(1,Math.min(5,Math.ceil(density*5)));
   for(let particle=0;particle<count;particle++){{const seed=(+frame.time_h*1000)+species*100000+index*17+particle*31;particles.push({{x:x+0.15+seeded(seed),y:y+0.15+seeded(seed+1),species}});}}
  }});
 }}
 return particles.slice(0,maxParticles);
}}
function drawAgents(frame){{
 const agents=frame.agents||[]; ctx.save(); ctx.globalAlpha=1;
 for(const agent of agents){{const x=agent.x+0.5,y=agent.y+0.5;ctx.fillStyle=agent.color||'#f8fafc';ctx.strokeStyle='#020617';ctx.lineWidth=.18;ctx.beginPath();ctx.arc(x,y,.42,0,Math.PI*2);ctx.fill();ctx.stroke();}}
 ctx.restore(); return agents;
}}
function drawParticles(frame){{
 const particles=populationParticles(frame); ctx.save(); ctx.globalAlpha=0.9;
 for(const particle of particles){{ctx.fillStyle=colors[particle.species%colors.length];ctx.beginPath();ctx.arc(particle.x,particle.y,.28,0,Math.PI*2);ctx.fill();}}
 const agents=drawAgents(frame); ctx.strokeStyle='#f8fafc'; ctx.lineWidth=.12;
 for(const agent of agents){{const previous=data.frames[Math.max(0,+slider.value-1)].agents||[];const match=previous.find(item=>item.id===agent.id);if(match){{ctx.beginPath();ctx.moveTo(match.x+.5,match.y+.5);ctx.lineTo(agent.x+.5,agent.y+.5);ctx.stroke()}}}}
 ctx.restore();
}}
function draw(){{
 const i=+slider.value, frame=data.frames[i], layer=selectedLayer(), field=fieldFor(frame,layer), flat=field.flat(), max=(layer.kind==='limit'||layer.kind==='phenotype')?1:Math.max(...flat.map(Number),1e-12);
 canvas.style.width=`${{data.width*+zoom.value}}px`; canvas.style.height=`${{data.height*+zoom.value}}px`; canvas.style.imageRendering=renderMode.value==='cells'?'pixelated':'auto'; if(renderMode.value!=='agents') drawField(field,layer,max); else {{ctx.clearRect(0,0,data.width,data.height);ctx.fillStyle='#080b10';ctx.fillRect(0,0,data.width,data.height)}} if(renderMode.value==='combined'||renderMode.value==='agents') drawParticles(frame);
 label.textContent=`t=${{frame.time_h.toFixed(2)}} h · frame ${{i+1}}/${{data.frames.length}}`; status.textContent=playing?'playing':'paused'; legend.textContent=(layer.kind==='limit'||layer.kind==='phenotype')?'categorical':`range 0-${{max.toPrecision(4)}}`;
 stats.innerHTML=Object.entries(frame.statistics.species_total_biomass).map(([name,value],j)=>`<span class="metric"><span class="swatch" style="background:${{colors[j%colors.length]}}"></span>${{name}} ${{value.toFixed(3)}} · ${{frame.statistics.species_occupied_cells[name]}} cells</span>`).join('')+`<span class="metric">mutations ${{frame.statistics.mutation_count}}</span>`;
 for(let s=0;s<frame.mutations.length;s++)for(let y=0;y<data.height;y++)for(let x=0;x<data.width;x++)if(frame.mutations[s][y][x]){{ctx.strokeStyle='#f4c95d';ctx.strokeRect(x+.15,y+.15,.7,.7);}}
 if(layer.kind==='rate') for(let y=0;y<data.height;y++)for(let x=0;x<data.width;x++) if(field[y][x]>0){{ctx.strokeStyle='#22c55e';ctx.strokeRect(x+.15,y+.15,.7,.7);}}
}}
canvas.onclick=(event)=>{{const rect=canvas.getBoundingClientRect(),x=Math.min(data.width-1,Math.floor((event.clientX-rect.left)*data.width/rect.width)),y=Math.min(data.height-1,Math.floor((event.clientY-rect.top)*data.height/rect.height)),frame=data.frames[+slider.value],layer=selectedLayer(), conditions=Object.entries(frame.conditions||{{}}).map(([name,values])=>`${{name}}=${{Number(values[y][x]).toPrecision(4)}}`).join(', ')||'no extra conditions'; let response=''; if(layer.kind==='rate'||layer.kind==='limit') response=` · rate=${{frame.niche.effective_growth_rate_per_h[names[layer.index]][y][x].toPrecision(4)}}/h · limiting=${{frame.niche.limiting_factor[names[layer.index]][y][x]}}`; if(layer.kind==='phenotype') response=` · state=${{frame.phenotypes[layer.index][y][x]}}`; inspect.textContent=`cell (${{x}},${{y}}): ${{conditions}}${{response}}`;}}
function setFrame(value){{slider.value=Math.max(0,Math.min(data.frames.length-1,value));draw()}} function toggle(){{playing=!playing;document.getElementById('play').textContent=playing?'⏸ Pause':'▶ Play';if(playing) timer=setInterval(()=>setFrame(+slider.value+1>=data.frames.length?0:+slider.value+1),+speed.value);else clearInterval(timer);draw()}}
slider.oninput=draw; selector.onchange=draw; renderMode.onchange=draw; zoom.oninput=draw; opacity.oninput=draw; speed.onchange=()=>{{if(playing){{clearInterval(timer);timer=setInterval(()=>setFrame(+slider.value+1>=data.frames.length?0:+slider.value+1),+speed.value)}}}}; document.getElementById('play').onclick=toggle; document.getElementById('back').onclick=()=>setFrame(+slider.value-1); document.getElementById('forward').onclick=()=>setFrame(+slider.value+1); document.addEventListener('keydown',event=>{{if(event.target.tagName==='INPUT'||event.target.tagName==='SELECT')return;if(event.code==='Space'){{event.preventDefault();toggle()}}if(event.key==='ArrowLeft')setFrame(+slider.value-1);if(event.key==='ArrowRight')setFrame(+slider.value+1);if(event.key==='Home')setFrame(0);if(event.key==='End')setFrame(data.frames.length-1)}}); draw();
</script>"""
    viewer_path.write_text(html, encoding="utf-8")
    return viewer_path
