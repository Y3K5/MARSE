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
canvas{{image-rendering:pixelated;border:1px solid #526071;max-width:80vw}}
button,input{{margin:4px}} label{{margin-right:12px}} .mut{{color:#f4c95d}}
</style>
<h1>MARSE ecosystem: {result.config.experiment_id}</h1>
<div><button id="play">Play</button><input id="time" type="range" min="0" max="{len(result.frames) - 1}" value="0" style="width:50%">
<span id="label"></span></div><div id="layers"></div><canvas id="view" width="{result.config.width}" height="{result.config.height}"></canvas>
<p class="mut">Color intensity shows the selected layer; mutations are outlined in yellow.</p>
<script>
const data = {data.read_text(encoding="utf-8")};
const canvas=document.getElementById('view'), ctx=canvas.getContext('2d');
const slider=document.getElementById('time'), label=document.getElementById('label'), layers=document.getElementById('layers');
let playing=false, timer;
const colors=['#46b3ff','#ff6b6b','#7ee081','#c084fc','#fb923c','#f472b6'];
const names=[...data.species], nutrientNames=[...data.nutrients];
for(const [kind, values] of [['species',names],['nutrient',nutrientNames]]) for(const [i,name] of values.entries()){{
 const id=kind+i; layers.innerHTML+=`<label><input type="radio" name="layer" id="${{id}}" ${{i===0?'checked':''}}>${{name}} (${{kind}})</label>`;
}}
function draw(){{
 const i=+slider.value, frame=data.frames[i], selected=document.querySelector('input[name=layer]:checked').id;
 const kind=selected.startsWith('species')?'species':'nutrients', index=+selected.match(/\\d+$/)[0], field=frame[kind][(kind==='species'?names:nutrientNames)[index]];
 const image=ctx.createImageData(data.width,data.height); let max=Math.max(...field.flat(),1e-12);
 for(let y=0;y<data.height;y++)for(let x=0;x<data.width;x++){{const p=(y*data.width+x)*4,v=Math.min(1,field[y][x]/max), c=kind==='species'?hex(colors[index%colors.length]):[80,220,150]; if(typeof c==='string') c=c.slice(1).match(/.{{2}}/g).map(x=>parseInt(x,16)); image.data.set([...c.map(x=>x*v),255],p);}}
 ctx.putImageData(image,0,0); label.textContent=`t=${{frame.time_h.toFixed(2)}} h`;
 for(let s=0;s<frame.mutations.length;s++)for(let y=0;y<data.height;y++)for(let x=0;x<data.width;x++)if(frame.mutations[s][y][x]){{ctx.strokeStyle='#f4c95d';ctx.strokeRect(x+.15,y+.15,.7,.7);}}
}}
function hex(s){{return s.match(/.{{2}}/g).map(x=>parseInt(x,16));}}
slider.oninput=draw; layers.onchange=draw; document.getElementById('play').onclick=()=>{{playing=!playing;document.getElementById('play').textContent=playing?'Pause':'Play'; if(playing) timer=setInterval(()=>{{slider.value=(+slider.value+1)%data.frames.length;draw()}},100);else clearInterval(timer)}}; draw();
</script>"""
    viewer_path.write_text(html, encoding="utf-8")
    return viewer_path
