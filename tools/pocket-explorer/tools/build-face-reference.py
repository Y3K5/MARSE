"""Derive a CC0 MakeHuman head locally. No source program code is imported.
Inputs are pinned public graphical assets; all registration choices are authored.
"""

import hashlib
import json
from pathlib import Path

import bpy

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "assets/face-reference"
vertices, faces, groups = [], [], []
group = ""
for line in (ASSETS / "base.obj").read_text().splitlines():
    if line.startswith("v "):
        vertices.append(list(map(float, line.split()[1:4])))
    elif line.startswith("g "):
        group = line[2:]
    elif line.startswith("f ") and group == "body":
        faces.append([int(s.split("/")[0]) - 1 for s in line.split()[1:]])
for name in ["african", "asian", "caucasian"]:
    for line in (ASSETS / (name + "-mouth-open.target")).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        index, *delta = line.split()
        for axis in range(3):
            vertices[int(index)][axis] += float(delta[axis]) * 1.3 / 3
# Only graphical head/neck surfaces enter the viewer; no body or helpers.
faces = [f for f in faces if all(vertices[i][1] > 5.84 for i in f)]
used = sorted({i for face in faces for i in face})
lookup = {v: i for i, v in enumerate(used)}
positions = [
    (vertices[i][0] * 80, (vertices[i][1] - 6.58) * 90 + 27, (vertices[i][2] - 1.50) * 80 + 46)
    for i in used
]
mesh = bpy.data.meshes.new("MakeHuman CC0 head")
mesh.from_pydata(positions, [], [[lookup[i] for i in f] for f in faces])
mesh.update()
obj = bpy.data.objects.new("Generic adult facial context", mesh)
bpy.context.collection.objects.link(obj)
bpy.context.view_layer.objects.active = obj
obj.select_set(True)
sub = obj.modifiers.new("Facial surface subdivision", "SUBSURF")
sub.levels = 2
bpy.ops.object.modifier_apply(modifier=sub.name)
obj.data.calc_loop_triangles()
packed = {
    "positions": [round(v, 5) for p in obj.data.vertices for v in p.co],
    "indices": [i for t in obj.data.loop_triangles for i in t.vertices],
}
inputs = {
    p.name: hashlib.sha256(p.read_bytes()).hexdigest()
    for p in ASSETS.iterdir()
    if p.suffix in [".obj", ".target"] or p.name == "LICENSE.ASSETS.md"
}
record = {
    "format": "marse.facial-context/1",
    "license": "CC0-1.0",
    "upstream": "https://github.com/makehumancommunity/makehuman",
    "commit": "a8bc2d54ff0ac92e78ff71431b1023eda42bf482",
    "inputs": inputs,
    "triangles": len(packed["indices"]) // 3,
    "pose": "Equal blend of three stock mouth-open targets, multiplier 1.3; artistic pose",
    "registration": {
        "scale": [80, 90, 80],
        "source_landmark": [0, 6.58, 1.50],
        "target_landmark": [0, 27, 46],
        "status": "authored, not anatomically measured",
    },
    "claim": "Graphical face asset; no clinical anatomy or patient identity claim",
    **packed,
}
(ROOT / "src/pocket-face-input.json").write_text(json.dumps(record, separators=(",", ":")) + "\n")
print(
    json.dumps({k: record[k] for k in ["format", "license", "commit", "triangles", "registration"]})
)
