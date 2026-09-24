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
body{{font:15px system-ui;background:#101318;color:#e7edf5;margin:24px}}
canvas{{image-rendering:pixelated;border:1px solid #526071;max-width:80vw;cursor:crosshair}}
button,input,select{{margin:4px}} label{{margin-right:12px}} .mut{{color:#f4c95d}}
#controls{{display:flex;flex-wrap:wrap;align-items:center;gap:6px}} #legend{{font-family:monospace}}
#inspect{{min-height:1.5em;color:#b8c7d9}}
</style>
<h1>MARSE ecosystem: {result.config.experiment_id}</h1>
<div id="controls"><button id="play">Play</button><input id="time" type="range" min="0" max="{len(result.frames) - 1}" value="0" style="width:50%">
<span id="label"></span><label>Layer <select id="layer"></select></label><span id="legend"></span></div>
<div id="stats"></div><canvas id="view" width="{result.config.width}" height="{result.config.height}"></canvas>
<div id="inspect">Click a cell to inspect its local conditions and niche response.</div>
<p class="mut">Niche-rate layers show effective growth; limiting-factor layers show the strongest constraint. Viable cells are outlined in green.</p>
<script>
const data = {data.read_text(encoding="utf-8")};
const canvas=document.getElementById('view'), ctx=canvas.getContext('2d'), stats=document.getElementById('stats'), inspect=document.getElementById('inspect');
const slider=document.getElementById('time'), label=document.getElementById('label'), selector=document.getElementById('layer'), legend=document.getElementById('legend');
let playing=false, timer;
const colors=['#46b3ff','#ff6b6b','#7ee081','#c084fc','#fb923c','#f472b6'];
const names=[...data.species], nutrientNames=[...data.nutrients], conditionNames=[...(data.conditions||[])];
const layers=[];
for(const [kind, values] of [['species',names],['nutrient',nutrientNames],['condition',conditionNames]]) for(const [i,name] of values.entries()) layers.push({{kind,index:i,label:`${{name}} (${{kind}})`}});
for(const [i,name] of names.entries()) layers.push({{kind:'rate',index:i,label:`${{name}} (effective growth rate)`}});
for(const [i,name] of names.entries()) layers.push({{kind:'limit',index:i,label:`${{name}} (limiting factor)`}});
for(const [i,layer] of layers.entries()) selector.innerHTML+=`<option value="${{i}}">${{layer.label}}</option>`;
function selectedLayer(){{return layers[+selector.value]}}
function fieldFor(frame, layer){{
 if(layer.kind==='species') return frame.species[names[layer.index]];
 if(layer.kind==='nutrient') return frame.nutrients[nutrientNames[layer.index]];
 if(layer.kind==='condition') return frame.conditions[conditionNames[layer.index]];
 if(layer.kind==='rate') return frame.niche.effective_growth_rate_per_h[names[layer.index]];
 return frame.niche.limiting_factor[names[layer.index]];
}}
function color(value, layer, max){{
 if(layer.kind==='limit'){{const palette=['#ef4444','#f59e0b','#22c55e','#38bdf8','#c084fc']; return palette[Math.abs(value.split('').reduce((a,c)=>a+c.charCodeAt(0),0))%palette.length]}}
 const ratio=Math.max(0,Math.min(1,value/Math.max(max,1e-12))); return `rgb(${{Math.round(40+180*ratio)}},${{Math.round(90+130*ratio)}},${{Math.round(150+90*ratio)}})`;
}}
function draw(){{
 const i=+slider.value, frame=data.frames[i], layer=selectedLayer(), field=fieldFor(frame,layer), flat=field.flat(), max=layer.kind==='limit'?1:Math.max(...flat.map(Number),1e-12);
 ctx.clearRect(0,0,data.width,data.height);
 for(let y=0;y<data.height;y++)for(let x=0;x<data.width;x++){{ctx.fillStyle=color(field[y][x],layer,max);ctx.fillRect(x,y,1,1);}}
 label.textContent=`t=${{frame.time_h.toFixed(2)}} h`;  legend.textContent=layer.kind==='limit'?'categorical limiting factors':`range 0-${{max.toPrecision(4)}}`;
 stats.textContent=Object.entries(frame.statistics.species_total_biomass).map(([name,value])=>`${{name}} biomass=${{value.toFixed(3)}} occupied=${{frame.statistics.species_occupied_cells[name]}}`).join(' | ')+` | mutations=${{frame.statistics.mutation_count}}`;
 for(let s=0;s<frame.mutations.length;s++)for(let y=0;y<data.height;y++)for(let x=0;x<data.width;x++)if(frame.mutations[s][y][x]){{ctx.strokeStyle='#f4c95d';ctx.strokeRect(x+.15,y+.15,.7,.7);}}
 if(layer.kind==='rate') for(let y=0;y<data.height;y++)for(let x=0;x<data.width;x++) if(field[y][x]>0){{ctx.strokeStyle='#22c55e';ctx.strokeRect(x+.15,y+.15,.7,.7);}}
}}
canvas.onclick=(event)=>{{const rect=canvas.getBoundingClientRect(),x=Math.min(data.width-1,Math.floor((event.clientX-rect.left)*data.width/rect.width)),y=Math.min(data.height-1,Math.floor((event.clientY-rect.top)*data.height/rect.height)),frame=data.frames[+slider.value],layer=selectedLayer(), conditions=Object.entries(frame.conditions||{{}}).map(([name,values])=>`${{name}}=${{values[y][x].toPrecision(4)}}`).join(', ')||'no extra conditions'; let response=''; if(layer.kind==='rate'||layer.kind==='limit') response=` rate=${{frame.niche.effective_growth_rate_per_h[names[layer.index]][y][x].toPrecision(4)}}/h, limiting=${{frame.niche.limiting_factor[names[layer.index]][y][x]}}`; inspect.textContent=`cell (${{x}},${{y}}): ${{conditions}}${{response}}`;}}
slider.oninput=draw; selector.onchange=draw; document.getElementById('play').onclick=()=>{{playing=!playing;document.getElementById('play').textContent=playing?'Pause':'Play'; if(playing) timer=setInterval(()=>{{slider.value=(+slider.value+1)%data.frames.length;draw()}},100);else clearInterval(timer)}}; draw();
</script>"""
    viewer_path.write_text(html, encoding="utf-8")
    return viewer_path
