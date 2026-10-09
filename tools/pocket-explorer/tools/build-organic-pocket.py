"""Local Blender authoring pipeline. Generic educational geometry, not measured anatomy.
Run with: blender --background --factory-startup --python tools/build-organic-pocket.py
"""

import json
import math
from itertools import pairwise
from pathlib import Path

import bmesh
import bpy
from mathutils import Vector
from mathutils.bvhtree import BVHTree
from mathutils.geometry import tessellate_polygon

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets/pocket-organic"
OUT.mkdir(parents=True, exist_ok=True)
# Identity of the earlier illustrative scaffold, retained as reference only.
# This source package rebuilds all geometry without that bulk reference file.
REFERENCE_SCAFFOLD_SHA256 = "944c7aea8ae12b1504413afd9016af52ce799221a9d65b2f3c2e32d17a85a186"
# Preserve reference identity; live geometry below uses explicitly authored rules.
bpy.ops.object.select_all(action="SELECT")
bpy.ops.object.delete(use_global=False)


def coords(p):
    return (p[0], -p[2], p[1])


def make(name, tid, pts, faces, preserve_orientation=False):
    m = bpy.data.meshes.new(name)
    m.from_pydata([coords(p) for p in pts], [], faces)
    m.update()
    obj = bpy.data.objects.new(name, m)
    bpy.context.collection.objects.link(obj)
    obj["tissueId"] = tid
    bm = bmesh.new()
    bm.from_mesh(m)
    bmesh.ops.remove_doubles(bm, verts=list(bm.verts), dist=0.000003)
    if not preserve_orientation:
        bmesh.ops.recalc_face_normals(bm, faces=list(bm.faces))
    bm.to_mesh(m)
    bm.free()
    for p in m.polygons:
        p.use_smooth = True
    return obj


def from_input(item):
    pts = [item["positions"][i : i + 3] for i in range(0, len(item["positions"]), 3)]
    ind = item["indices"] or list(range(len(pts)))
    return make(
        item["name"],
        item["id"],
        pts,
        [ind[i : i + 3] for i in range(0, len(ind), 3)],
        item["id"] in ["enamel", "cementum", "pdl"],
    )


def apply(obj, mod):
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    bpy.ops.object.modifier_apply(modifier=mod.name)
    obj.select_set(False)


def boolean(obj, cutter, operation="DIFFERENCE"):
    bpy.context.view_layer.update()
    print("CSG", obj.name, operation, cutter.name, flush=True)
    m = obj.modifiers.new("Anatomical compartment boundary", "BOOLEAN")
    m.operation = operation
    m.solver = "MANIFOLD"
    m.object = cutter
    apply(obj, m)


def sphere(name, tid, c, r):
    bpy.ops.mesh.primitive_uv_sphere_add(segments=80, ring_count=48, location=coords(c))
    o = bpy.context.object
    o.name = name
    o["tissueId"] = tid
    o.scale = (r[0], r[2], r[1])
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    for p in o.data.polygons:
        p.use_smooth = True
    o.select_set(False)
    return o


def combine(objects, name, tid, voxel=None):
    bpy.ops.object.select_all(action="DESELECT")
    for o in objects:
        o.select_set(True)
    bpy.context.view_layer.objects.active = objects[0]
    bpy.ops.object.join()
    o = bpy.context.object
    o.name = name
    o["tissueId"] = tid
    if voxel:
        m = o.modifiers.new("Continuous dentin volume", "REMESH")
        m.mode = "VOXEL"
        m.voxel_size = voxel
        m.use_smooth_shade = True
        apply(o, m)
        smooth = o.modifiers.new("Organic transition", "SMOOTH")
        smooth.factor = 0.35
        smooth.iterations = 3
        apply(o, smooth)
        m = o.modifiers.new("Browser mesh budget", "DECIMATE")
        m.ratio = 0.40 if tid in ["gingiva", "connective"] else 0.52
        apply(o, m)
        # The voxel union can leave microscopic disconnected remeshing slivers.
        # Remove only small authoring artifacts; retain the complete main body.
        bm = bmesh.new()
        bm.from_mesh(o.data)
        unseen = set(bm.verts)
        groups = []
        while unseen:
            group = {unseen.pop()}
            todo = list(group)
            while todo:
                v = todo.pop()
                for edge in v.link_edges:
                    other = edge.other_vert(v)
                    if other in unseen:
                        unseen.remove(other)
                        group.add(other)
                        todo.append(other)
            groups.append(group)
        groups.sort(key=len, reverse=True)
        membership = {v: i for i, g in enumerate(groups) for v in g}
        volumes = [0.0] * len(groups)
        for f in bm.faces:
            volumes[membership[f.verts[0]]] += (
                f.calc_center_median().dot(f.normal) * f.calc_area() / 3
            )
        # These authored volumes have no intended enclosed micro-cavities.
        # Tiny detached positive scraps and negative remeshing bubbles are
        # repaired together; substantial secondary volumes remain fatal.
        discarded = list(range(1, len(groups)))
        minor = set().union(*(groups[i] for i in discarded)) if discarded else set()
        fraction = sum(abs(volumes[i]) for i in discarded) / max(abs(volumes[0]), 1e-12)
        print(
            "REMESH_COMPONENTS",
            name,
            "main",
            len(groups[0]),
            "positive_scrap_fraction",
            fraction,
            "cavities",
            sum(v < 0 for v in volumes),
            flush=True,
        )
        if tid == "bone":
            assert fraction < 0.001, "Substantial detached bone volume requires review"
        elif tid in ["gingiva", "connective"]:
            assert fraction < 0.005, "Substantial detached soft-tissue volume requires review"
        else:
            assert all(len(groups[i]) <= 12 for i in discarded), (
                "Substantial disconnected tooth requires review"
            )
        o["removed_remesh_scrap_volume_fraction"] = fraction
        o["filled_artificial_cavity_components"] = sum(volumes[i] < 0 for i in discarded)
        o["removed_remesh_sliver_components"] = len(discarded)
        o["removed_remesh_sliver_vertices"] = len(minor)
        if minor:
            bmesh.ops.delete(bm, geom=list(minor), context="VERTS")
        # Close only tiny remeshing pinholes; larger defects remain fatal.
        boundary = {e for e in bm.edges if e.is_boundary}
        filled_count = 0
        max_span = 0
        while boundary:
            component = {boundary.pop()}
            todo = list(component)
            while todo:
                e = todo.pop()
                for v in e.verts:
                    for q in v.link_edges:
                        if q in boundary:
                            boundary.remove(q)
                            component.add(q)
                            todo.append(q)
            vs = {v for e in component for v in e.verts}
            span = max(max(v.co[i] for v in vs) - min(v.co[i] for v in vs) for i in range(3))
            assert len(component) <= 12, "Large remeshing hole requires reauthoring"
            assert span <= voxel, "Large remeshing hole requires reauthoring"
            bmesh.ops.holes_fill(bm, edges=list(component), sides=12)
            filled_count += 1
            max_span = max(max_span, span)
        o["filled_remesh_pinholes"] = filled_count
        o["maximum_pinhole_span"] = max_span
        wire = [e for e in bm.edges if not e.link_faces]
        o["removed_remesh_wire_edges"] = len(wire)
        if wire:
            bmesh.ops.delete(bm, geom=wire, context="EDGES")
        isolated = [v for v in bm.verts if not v.link_edges]
        if isolated:
            bmesh.ops.delete(bm, geom=isolated, context="VERTS")
        bm.normal_update()
        print(
            "EDGE_TYPES",
            name,
            {n: sum(len(e.link_faces) == n for e in bm.edges) for n in [0, 1, 3, 4]},
            flush=True,
        )
        bm.to_mesh(o.data)
        bm.free()
    o.select_set(False)
    return o


def clone(obj, name, tid=None):
    bpy.context.view_layer.update()
    o = bpy.data.objects.new(name, obj.data.copy())
    bpy.context.collection.objects.link(o)
    o.matrix_world = obj.matrix_world.copy()
    o["tissueId"] = tid or obj["tissueId"]
    o["partName"] = name
    for k in obj.keys():  # noqa: SIM118 -- Blender ID properties are not iterable
        if k.startswith(("removed_", "filled_", "maximum_", "retained_")) or k == "repair":
            o[k] = obj[k]
    bpy.context.view_layer.update()
    return o


def discard(o):
    bpy.data.objects.remove(o, do_unlink=True)


def smoothstep(t):
    t = max(0, min(1, t))
    return t * t * (3 - 2 * t)


def sp(x, p):
    return math.copysign(abs(x) ** p, x)


# Mesiobuccal, distobuccal, distal, mesiolingual and distolingual cusps.
# Lingual cusps stand higher than buccal cusps, as on mandibular molars.
CUSPS = [
    (-2.55, 2.10, 1.30),
    (0.65, 2.30, 1.18),
    (3.30, 1.20, 0.84),
    (-2.30, -2.25, 1.58),
    (1.70, -2.10, 1.44),
]
# Mesiobuccal, distobuccal and lingual developmental grooves (offset, side).
GROOVES = [(-0.85, 1), (2.0, 1), (-0.10, -1)]
# Authored crown descriptions. Half widths are mesiodistal and buccolingual.
MOLAR = {
    "cusps": CUSPS,
    "grooves": GROOVES,
    "base": 4.30,
    "half": (5.10, 4.25),
    "taper": 0.07,
    "lingual_exponent": 0.70,
    "ridges": (-4.05, 3.95),
    "fossae": ((-3.25, 0.0), (2.55, 0.05)),
    "central_span": 17.0,
}
# Mandibular second premolar (three-cusp form): one dominant buccal cusp.
PREMOLAR = {
    "cusps": [(0.0, 1.30, 1.55), (-1.00, -1.55, 1.05), (1.05, -1.60, 0.85)],
    "grooves": [(0.02, -1)],
    "base": 4.55,
    "half": (3.55, 4.00),
    "taper": 0.03,
    "lingual_exponent": 0.78,
    "ridges": (-2.75, 2.75),
    "fossae": ((-1.85, -0.15), (1.85, -0.20)),
    "central_span": 4.0,
}


def occlusal(x, z, grooves=1.0, spec=MOLAR):
    h = spec["base"]
    for cx, cz, amp in spec["cusps"]:
        dx = x - cx
        dz = z - cz
        h += amp * math.exp(-(dx * dx / 2.45 + dz * dz / 2.15))
        # Triangular ridge running from each cusp tip toward the central groove.
        vx = -0.10 * cx
        vz = -cz
        length = math.hypot(vx, vz)
        u = (dx * vx + dz * vz) / length
        w = (dx * vz - dz * vx) / length
        along = smoothstep((u + 0.25) / 0.45) * math.exp(-((u / (0.75 * length)) ** 2))
        h += 0.10 * amp * math.exp(-w * w / 0.14) * along
    # Mesial and distal marginal ridges close the occlusal table.
    for mx in spec["ridges"]:
        h += 0.20 * math.exp(-(((x - mx) / 0.50) ** 2)) * math.exp(-((z / 2.6) ** 4))
    # Authored central and developmental fissures, not caries or scan data.
    central = z + 0.12 - 0.18 * math.sin(x * 0.9)
    h -= (
        grooves
        * 0.50
        * math.exp(-central * central / 0.045)
        * math.exp(-x * x / spec["central_span"])
    )
    for q, sign in spec["grooves"]:
        line = x - q - 0.15 * z
        reach = smoothstep((sign * z + 0.1) / 0.6)
        h -= grooves * 0.30 * math.exp(-line * line / 0.042) * reach
    # Mesial and distal triangular fossae at the ends of the central groove.
    for fx, fz in spec["fossae"]:
        h -= grooves * 0.22 * math.exp(-((x - fx) ** 2 + (z - fz) ** 2) / 0.28)
    return h


def cej_height(a):
    # Cervical line curves occlusally on the proximal surfaces.
    return -0.12 + 0.12 * math.cos(2 * a)


def crown_outline(a, spec=MOLAR):
    # Squarer mesial outline, distal taper, and a flatter lingual than buccal face.
    hx, hz = spec["half"]
    s = math.sin(a)
    x = hx * sp(math.cos(a), 0.78) * (1 - spec["taper"] * math.cos(a))
    z = hz * sp(s, 0.83 if s > 0 else spec["lingual_exponent"])
    return x, z


def axial_scale(a, t):
    s = math.sin(a)
    # Buccal height of contour in the cervical third, lingual in the middle third;
    # the buccal surface converges more strongly toward the occlusal table, and
    # the walls round over into the cusp slopes instead of meeting them at a rim.
    peak = 0.40 - 0.13 * s
    top = 0.830 - 0.055 * s
    if t < peak:
        return 0.862 + 0.138 * math.sin(0.5 * math.pi * t / peak)
    return 1 - (1 - top) * ((t - peak) / (1 - peak)) ** 2.4


def angle_gap(a, b):
    return math.atan2(math.sin(a - b), math.cos(a - b))


def crown_mesh(inner=False, spec=MOLAR, name=None):
    na = 144
    rows = []
    for j in range(49):
        t = j / 48
        # Height rises quickly through the cervical and middle thirds, then eases
        # so the axial wall turns smoothly into the outer cusp slopes.
        rise = t + 0.55 * (t - t * t)
        row = []
        for i in range(na):
            a = i / na * math.tau
            ox, oz = crown_outline(a, spec)
            k = axial_scale(a, t)
            # Each cusp swells its part of the occlusal third into a lobe; the
            # developmental grooves continue between lobes onto the axial faces.
            lobe = sum(
                0.035 * math.exp(-((angle_gap(a, math.atan2(cz, cx)) / 0.38) ** 2))
                for cx, cz, _ in spec["cusps"]
            )
            k *= 1 + lobe * smoothstep((t - 0.45) / 0.40)
            dent = 0
            for q, sign in spec["grooves"]:
                line = ox * k - q - 0.15 * oz * k
                face = smoothstep((sign * math.sin(a) - 0.35) / 0.3)
                dent += 0.07 * math.exp(-line * line / 0.07) * face
            k *= 1 - dent * smoothstep((t - 0.42) / 0.40) / max(math.hypot(ox, oz), 1)
            x = ox * k
            z = oz * k
            cej = cej_height(a)
            y = cej + rise * (occlusal(x, z, 1.0, spec) - cej)
            if inner:
                thickness = 0.014 + 0.14 * smoothstep(t / 0.6)
                x *= 1 - thickness
                z *= 1 - thickness
                y -= 0.045 + 0.62 * smoothstep(t)
            row.append((x, y, z))
        rows.append(row)
    rim = rows[-1]
    for k in range(1, 35):
        r = 1 - k / 35
        row = []
        for x, _, z in rim:
            xx = x * r
            zz = z * r
            yy = (
                occlusal(xx / 0.846, zz / 0.846, 0.3, spec) - 0.665
                if inner
                else occlusal(xx, zz, 1.0, spec)
            )
            row.append((xx, yy, zz))
        rows.append(row)
    pts = [p for r in rows for p in r]
    faces = []
    for j in range(len(rows) - 1):
        for i in range(na):
            faces.append(
                (j * na + i, (j + 1) * na + i, (j + 1) * na + (i + 1) % na, j * na + (i + 1) % na)
            )
    bot = len(pts)
    pts.append((0, -0.19 if not inner else -0.23, 0))
    tip = len(pts)
    pts.append((0, occlusal(0, 0, 0.3, spec) - 0.665 if inner else occlusal(0, 0, 1.0, spec), 0))
    for i in range(na):
        faces.extend(
            [
                (bot, i, (i + 1) % na),
                ((len(rows) - 1) * na + i, tip, (len(rows) - 1) * na + (i + 1) % na),
            ]
        )
    return make(
        name or ("Coronal dentin" if inner else "Five-cusp enamel envelope"),
        "dentin" if inner else "enamel",
        pts,
        faces,
    )


def root_center(side, t):
    return (
        (-2.05 - 1.0 * math.sin(t * math.pi / 2) + 1.45 * t * t)
        if side < 0
        else (2.0 + 1.05 * math.sin(t * math.pi / 2) + 0.85 * t * t),
        -2 - (12.45 if side < 0 else 11.60) * t,
        0.10 * math.sin(t * math.pi) * side,
    )


def root_taper(t):
    # Conical taper ending in a blunt, rounded apex rather than a needle point.
    t = max(0.0, min(1.0, t))
    return (1 - 0.50 * t) * (1 - t**3.2) ** 0.42


def root_radii(side, t):
    taper = root_taper(t)
    return (
        (1.62 if side < 0 else 1.74) * taper + 0.03,
        (2.15 if side < 0 else 1.80) * taper + 0.03,
    )


def rings_volume(name, tid, rows, na=64, concavity=None):
    pts = []
    faces = []
    for j, (x, y, z, rx, rz) in enumerate(rows):
        for i in range(na):
            a = i / na * math.tau
            dent = concavity(a, j / (len(rows) - 1)) if concavity else 0
            pts.append((x + rx * (1 - dent) * math.cos(a), y, z + rz * math.sin(a)))
    for j in range(len(rows) - 1):
        for i in range(na):
            faces.append(
                (j * na + i, (j + 1) * na + i, (j + 1) * na + (i + 1) % na, j * na + (i + 1) % na)
            )
    for end in [0, len(rows) - 1]:
        index = len(pts)
        pts.append(rows[end][:3])
        base = end * na
        for i in range(na):
            faces.append((index, base + i, base + (i + 1) % na))
    return make(name, tid, pts, faces)


def root(side):
    rows = [(*root_center(side, j / 75), *root_radii(side, j / 75)) for j in range(76)]
    # Proximal root concavities: both faces of the broad mesial root and the
    # furcal face of the distal root. The pocket-facing distal face stays convex.
    depth = {(-1, 1): 0.20, (-1, -1): 0.13, (1, -1): 0.11}

    def concavity(a, t):
        g = depth.get((side, 1 if math.cos(a) > 0 else -1), 0)
        along = smoothstep(t / 0.12) * (1 - smoothstep((t - 0.55) / 0.4))
        return g * math.exp(-(math.sin(a) ** 2) / 0.22) * along

    return rings_volume("Mesial root" if side < 0 else "Distal root", "dentin", rows, 72, concavity)


def shell_from(inner, outer, name, tid):
    inner.data.calc_loop_triangles()
    outer.data.calc_loop_triangles()
    pts = []
    faces = []
    for o, reverse in [(outer, False), (inner, True)]:
        offset = len(pts)
        pts.extend((q.x, q.z, -q.y) for q in [o.matrix_world @ v.co for v in o.data.vertices])
        faces.extend(
            tuple(offset + i for i in (reversed(t.vertices) if reverse else t.vertices))
            for t in o.data.loop_triangles
        )
    return make(name, tid, pts, faces, True)


def expanded(obj, distance, name, tid):
    o = clone(obj, name, tid)
    o.data.update()
    for v in o.data.vertices:
        v.co += v.normal * distance
    o.data.update()
    return o


def superellipsoid(name, tid, c, r, e=0.5, rings=36, segs=72):
    # Rounded box: squarer than an ellipsoid, without sharp edges.
    pts = []
    faces = []
    for j in range(1, rings):
        v = -math.pi / 2 + j / rings * math.pi
        cv = sp(math.cos(v), e)
        sv = sp(math.sin(v), e)
        for i in range(segs):
            u = i / segs * math.tau
            pts.append(
                (
                    c[0] + r[0] * cv * sp(math.cos(u), e),
                    c[1] + r[1] * sv,
                    c[2] + r[2] * cv * sp(math.sin(u), e),
                )
            )
    for j in range(rings - 2):
        for i in range(segs):
            faces.append(
                (
                    j * segs + i,
                    j * segs + (i + 1) % segs,
                    (j + 1) * segs + (i + 1) % segs,
                    (j + 1) * segs + i,
                )
            )
    bottom = len(pts)
    pts.append((c[0], c[1] - r[1], c[2]))
    top = len(pts)
    pts.append((c[0], c[1] + r[1], c[2]))
    last = (rings - 2) * segs
    for i in range(segs):
        faces.append((bottom, (i + 1) % segs, i))
        faces.append((top, last + i, last + (i + 1) % segs))
    return make(name, tid, pts, faces)


def pulp_body():
    # Chamber follows the crown outline: wider mesiodistally, roof near the
    # cervical third, floor close to the cervical line.
    objs = [
        superellipsoid(
            "Pulp chamber roof and floor", "pulp", (0.05, 0.95, -0.05), (2.45, 0.95, 1.80)
        )
    ]
    for cx, cz, amp in CUSPS:
        rows = []
        for j in range(30):
            t = j / 29
            rows.append(
                (
                    cx * (0.38 + 0.22 * t),
                    1.45 + t * (1.25 + amp * 0.28),
                    cz * (0.36 + 0.22 * t),
                    0.60 * (1 - t) ** 0.70 + 0.05,
                    0.48 * (1 - t) ** 0.70 + 0.045,
                )
            )
        objs.append(rings_volume("Pulp horn", "pulp", rows, 32))
    for side, shift in [(-1, -0.55), (-1, 0.55), (1, 0)]:
        rows = []
        for j in range(85):
            u = j / 84
            if u < 0.22:
                t = u / 0.22
                x = side * 0.6 + (side * 2.03 - side * 0.6) * smoothstep(t)
                y = 0.55 - 2.55 * t
                z = shift * smoothstep(t)
                r = 0.49 - 0.13 * t
            else:
                t = (u - 0.22) / 0.78
                x, y, z = root_center(side, t)
                z += shift * (1 - 0.76 * t)
                r = 0.36 * (1 - t) ** 0.9 + 0.040
            rows.append((x, y, z, r, r * (1.18 if side > 0 else 0.86)))
        objs.append(rings_volume("Curving tapered root canal", "pulp", rows, 32))
    return combine(objs, "Continuous pulp chamber horns and tapered canals", "pulp", 0.065)


# Three-tooth segment: FDI 35 (mesial, x < 0), the selected FDI 36 at the origin,
# and FDI 37 (distal). Neighbours are context; 37 is the 36 model scaled to 0.94.
# The browser trims every part at SEGMENT_ENDS, so cut faces show tissue layers.
BEND = 0.006
NEIGHBOURS = {"35": (-8.93, 1.0, PREMOLAR), "37": (9.92, 0.94, MOLAR)}
SEGMENT_ENDS = (-9.5, 10.4)
TEETH = [(-8.93, False), (0.0, True), (9.92, True)]


def bend(x):
    # Gentle arch curvature along the segment.
    return BEND * x * x


def crest(x, z, health):
    """Alveolar crest height: interdental septa higher than mid-buccal/lingual
    crests around each tooth; periodontitis loss stays at the 36 distal site."""
    num = 0.0
    den = 0.0
    lift = 0.0
    for cx, multi in TEETH:
        dx = x - cx
        dz = z - bend(cx)
        w = math.exp(-dx * dx / 12.25)
        num += w * 0.30 * (dx * dx - dz * dz) / (dx * dx + dz * dz + 0.5)
        den += w
        if multi:
            # Interradicular bone reaches the furcation entrance.
            lift += 0.35 * math.exp(-(dx * dx + dz * dz) / 4)
    y = -3.65 + num / den + lift
    if health == "periodontitis":
        r = math.hypot(x, z)
        a = math.atan2(z, x)
        y -= (
            5.15
            * (0.5 + 0.5 * math.cos(a)) ** 4
            * smoothstep((r - 2.5) / 1.5)
            * math.exp(-(max(0.0, r - 4.6) ** 2) / 4.84)
        )
    return y


def half_width(y, inset=0.0):
    # Alveolar process half-width, flaring slightly below the crest.
    t = (y + 17) / 13.35
    return (5.5 - inset) * (0.91 + 0.09 * math.sin(max(0.0, min(1.25, t)) * math.pi * 0.8))


def bone_section(x, health, inset=0.0):
    """Closed (z, y) cross-section of the alveolar process at one x."""
    b = bend(x)
    bottom = -17 + inset
    zb = b + half_width(-3.65, inset)
    zl = b - half_width(-3.65, inset)
    cb = crest(x, zb, health) - 0.8 * inset
    cl = crest(x, zl, health) - 0.8 * inset
    pts = [(z, bottom) for z in linspace(zl, zb, 10)[:-1]]
    pts += [(b + half_width(y, inset), y) for y in linspace(bottom, cb, 30)[:-1]]
    zb_top = b + half_width(cb, inset)
    zl_top = b - half_width(cl, inset)
    pts += [(z, crest(x, z, health) - 0.8 * inset) for z in linspace(zb_top, zl_top, 70)[:-1]]
    pts += [(b - half_width(y, inset), y) for y in linspace(cl, bottom, 30)[:-1]]
    return pts


def loft(name, tid, sections, xs):
    """Closed volume lofted through equal-length (z, y) sections at each x."""
    n = len(sections[0])
    pts = [(x, y, z) for x, sec in zip(xs, sections, strict=True) for z, y in sec]
    faces = []
    for i in range(len(xs) - 1):
        for j in range(n):
            faces.append(
                (i * n + j, i * n + (j + 1) % n, (i + 1) * n + (j + 1) % n, (i + 1) * n + j)
            )
    for i in (0, len(xs) - 1):
        poly = [Vector((z, y, 0)) for z, y in sections[i]]
        for tri in tessellate_polygon([poly]):
            faces.append(tuple(i * n + k for k in tri))
    return make(name, tid, pts, faces)


def bone_tube(health, inset=0.0):
    margin = 2.2 + (1.0 if inset else 0.0)
    xs = linspace(SEGMENT_ENDS[0] - margin, SEGMENT_ENDS[1] + margin, 210)
    return loft(
        "Trabecular alveolar interior" if inset else "Buccal and lingual cortical plates",
        "bone",
        [bone_section(x, health, inset) for x in xs],
        xs,
    )


def neighbour_cone(x, z):
    """Gingival surface rising toward each neighbour's scalloped margin."""
    y = -99.0
    for cx, s, spec in NEIGHBOURS.values():
        dx = (x - cx) / s
        dz = (z - bend(cx)) / s
        a = math.atan2(dz, dx)
        ox, oz = crown_outline(a, spec)
        d = max(0.0, (math.hypot(dx, dz) - 0.862 * math.hypot(ox, oz)) * s)
        y = max(y, (0.35 - 0.30 * math.sin(a) ** 2) * s - 0.55 * d - 0.45 * d * d)
    return y


def gum_top(x, z, health):
    # Attached gingiva over the crest, rising into the neighbours' margins; the
    # cones of two adjacent teeth meet in an interdental papilla.
    return smax(crest(x, z, health) + COVER_TOP, neighbour_cone(x, z), 0.4)


def slab_section(x, health, inset=0.0):
    """(z, y) section of the continuous gingiva across the alveolar process."""
    b = bend(x)
    deep = 0.15 if inset else 0.0
    zb = b + half_width(-3.65)
    zl = b - half_width(-3.65)
    cb = crest(x, zb, health)
    cl = crest(x, zl, health)
    zb_top = b + half_width(cb)
    zl_top = b - half_width(cl)
    gb = gum_top(x, zb_top, health)
    gl = gum_top(x, zl_top, health)
    base_b = smin(-7.3, cb - 2.0, 1.0)
    base_l = smin(-7.3, cl - 2.0, 1.0)

    def cover(y, top, base):
        return COVER_SIDE - (COVER_SIDE - COVER_BASE) * smoothstep((top - y) / (top - base))

    def roll(face, base, side):
        # Rolled basal edge entering the cortical plate steeply.
        wb = b + side * half_width(base)
        start = (face[0], face[1])
        return bezier(
            (wb + side * 0.30 - side * 0.60, base - 0.50),
            (wb + side * 0.25, base - 0.45),
            (start[0], base - 0.15),
            start,
            14,
        )

    buccal = [
        (b + half_width(y) + cover(y, cb, base_b), y) for y in linspace(base_b + 0.3, cb - 0.05, 30)
    ]
    corner_b = [
        (zb_top + COVER_SIDE * math.cos(t), cb + (gb - cb) * math.sin(t))
        for t in linspace(0, math.pi / 2, 12)
    ]
    top = [(z, gum_top(x, z, health)) for z in linspace(zb_top, zl_top, 110)]
    corner_l = [
        (zl_top - COVER_SIDE * math.cos(t), cl + (gl - cl) * math.sin(t))
        for t in linspace(math.pi / 2, 0, 12)
    ]
    lingual = [
        (b - half_width(y) - cover(y, cl, base_l), y) for y in linspace(cl - 0.05, base_l + 0.3, 30)
    ]
    free = (
        roll(buccal[0], base_b, 1)[:-1]
        + buccal
        + corner_b[1:]
        + top[1:-1]
        + corner_l[:-1]
        + lingual
        + roll(lingual[-1], base_l, -1)[::-1][1:]
    )
    free = resample(free, 260)
    if inset:
        free = offset_chain(free, -inset)
    end_l = free[-1]
    end_b = free[0]
    inner_l = [(b - half_width(y) + 0.6 + deep, y) for y in linspace(end_l[1], cl - 0.7 - deep, 12)]
    inner_top = [
        (z, crest(x, z, health) - 0.7 - deep)
        for z in linspace(zl_top + 0.6 + deep, zb_top - 0.6 - deep, 40)
    ]
    inner_b = [(b + half_width(y) - 0.6 - deep, y) for y in linspace(cb - 0.7 - deep, end_b[1], 12)]
    embedded = resample(inner_l + inner_top + inner_b, 60)
    section = free + embedded[1:-1]
    return section


def slab_volumes(health):
    xs = linspace(SEGMENT_ENDS[0] - 1.2, SEGMENT_ENDS[1] + 1.2, 230)
    out = []
    for inset, name, tid in [
        (0.0, "Continuous gingiva over the three-tooth segment", "gingiva"),
        (EPITHELIUM, "Subepithelial connective tissue", "connective"),
    ]:
        sections = [slab_section(x, health, inset) for x in xs]
        for k in (0, 57, 115, 172, 229):
            assert simple_polygon(sections[k]), ("Self-intersecting gum slab", health, xs[k], inset)
        out.append(loft(name, tid, sections, xs))
    return out


def surface_tree(objects):
    bpy.context.view_layer.update()
    verts = []
    polys = []
    for o in objects:
        offset = len(verts)
        verts.extend(o.matrix_world @ v.co for v in o.data.vertices)
        polys.extend([offset + i for i in p.vertices] for p in o.data.polygons)
    return BVHTree.FromPolygons(verts, polys)


def cast(tree, origin, direction):
    hit = tree.ray_cast(Vector(coords(origin)), Vector(coords(direction)).normalized())
    return None if hit[0] is None else hit[3]


def outer_radius(tree, a, y, far=26.0):
    # Outermost surface met by a horizontal ray travelling toward the tooth axis.
    ca, sa = math.cos(a), math.sin(a)
    dist = cast(tree, (far * ca, y, far * sa), (-ca, 0, -sa))
    return None if dist is None else far - dist


def smax(a, b, k):
    return max(a, b) + max(k - abs(a - b), 0) ** 2 / (4 * k)


def smin(a, b, k):
    return -smax(-a, -b, k)


def linspace(a, b, n):
    return [a + (b - a) * k / (n - 1) for k in range(n)]


def bezier(p0, p1, p2, p3, n=24):
    out = []
    for k in range(n):
        t = k / (n - 1)
        u = 1 - t
        out.append(
            tuple(
                u**3 * p0[i] + 3 * u * u * t * p1[i] + 3 * u * t * t * p2[i] + t**3 * p3[i]
                for i in range(2)
            )
        )
    return out


def resample(points, n):
    # Fixed-count arc-length resampling keeps vertex rows aligned between azimuths.
    lengths = [0.0]
    for p, q in pairwise(points):
        lengths.append(lengths[-1] + math.dist(p, q))
    out = []
    j = 0
    for k in range(n):
        s = lengths[-1] * k / (n - 1)
        while j < len(points) - 2 and lengths[j + 1] < s:
            j += 1
        span = lengths[j + 1] - lengths[j]
        f = 0 if span == 0 else (s - lengths[j]) / span
        out.append(tuple(points[j][i] + f * (points[j + 1][i] - points[j][i]) for i in range(2)))
    return out


def offset_chain(points, distance):
    # Move a smooth open (radius, height) chain toward the tissue side.
    out = []
    for i, p in enumerate(points):
        a = points[max(i - 1, 0)]
        b = points[min(i + 1, len(points) - 1)]
        tx, ty = b[0] - a[0], b[1] - a[1]
        length = math.hypot(tx, ty) or 1
        out.append((p[0] + distance * ty / length, p[1] - distance * tx / length))
    return out


def crossing(p, q, r, s):
    d = (q[0] - p[0]) * (s[1] - r[1]) - (q[1] - p[1]) * (s[0] - r[0])
    if abs(d) < 1e-12:
        return None
    t = ((r[0] - p[0]) * (s[1] - r[1]) - (r[1] - p[1]) * (s[0] - r[0])) / d
    u = ((r[0] - p[0]) * (q[1] - p[1]) - (r[1] - p[1]) * (q[0] - p[0])) / d
    return t if 0 <= t <= 1 and 0 <= u <= 1 else None


def simple_polygon(points):
    n = len(points)
    for i in range(n):
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            if crossing(points[i], points[(i + 1) % n], points[j], points[(j + 1) % n]) is not None:
                return False
    return True


# Illustrative epithelial thickness and soft-tissue cover, exaggerated for reading.
EPITHELIUM = 0.20
COVER_TOP = 0.78
COVER_SIDE = 0.52
COVER_BASE = 0.40
FURCATION_FORNIX = -3.55


def soft_landmarks(health, a):
    d = (0.5 + 0.5 * math.cos(a)) ** 4 if health == "periodontitis" else 0
    # Scalloped margin: highest at the proximal surfaces, lower mid-buccally and
    # mid-lingually. The selected distal site keeps the preset landmarks exactly.
    margin = 0.35 - 0.30 * math.sin(a) ** 2 - 0.50 * d
    attach = -1.65 - 4.55 * d
    return d, margin, attach


COLLAR_REACH = 3.0
HIDDEN = 0.3


def soft_sections(health, a, tooth):
    """Closed (radius, height) sections of the 36 gingival collar and its core.

    The free gingiva is authored against the ray-cast tooth surface and rises out
    of the continuous gingiva. Everything beyond the free-gingival slope runs
    HIDDEN below the continuous gingiva or inside hard tissue, so the union and
    the later socket, cementum and pocket cuts cross steeply, never grazing.
    """
    ca, sa = math.cos(a), math.sin(a)
    d, margin, attach = soft_landmarks(health, a)

    def at(r):
        return r * ca, r * sa

    def bone_top(r):
        return crest(*at(r), health)

    # Radius where the ray leaves the crest top through the buccal or lingual face.
    edge = 99.0
    if abs(sa) > 1e-3:
        lo, hi = 0.0, 40.0
        for _ in range(40):
            mid = (lo + hi) / 2
            x, z = at(mid)
            if abs(z - bend(x)) < half_width(-3.65):
                lo = mid
            else:
                hi = mid
        edge = lo
    grid = linspace(margin + 0.4, -11.5, 260)
    tooth_r = []
    for y in grid:
        r = outer_radius(tooth, a, y)
        tooth_r.append(r if r is not None and r > 1.2 else tooth_r[-1] if tooth_r else 4.6)

    def rt(y):
        f = (grid[0] - y) / (grid[0] - grid[-1]) * (len(grid) - 1)
        i = max(0, min(len(grid) - 2, int(f)))
        f -= i
        return tooth_r[i] + (tooth_r[i + 1] - tooth_r[i]) * f

    def gap(y):
        u = max(0.0, min(1.0, (y - attach) / (margin - attach)))
        return 0.11 + (0.03 + 0.14 * d) * math.sin(math.pi * u) ** 0.8

    wall_top = margin - 0.18
    rw = rt(wall_top) + gap(wall_top)
    p3 = (rw + 0.48, margin - 0.16)
    reach = min(rt(attach) + COLLAR_REACH, edge)
    # The surface the cone dives into: continuous gingiva, HIDDEN below its top,
    # then (toward the cheek or tongue) round its edge and down the plate.
    # Where the neighbouring gingiva stands higher than this margin (a receded
    # site beside a healthy papilla) the collar stays below its own margin and
    # the continuous gingiva forms the visible surface there.
    cover = [
        (r, min(gum_top(*at(r), health) - HIDDEN, p3[1] - 0.25 - 0.5 * (r - p3[0])))
        for r in linspace(p3[0], reach, 44)
    ]
    if reach >= edge:
        x, z = at(edge)
        top = crest(x, z, health)
        lift = cover[-1][1] - top
        side = COVER_SIDE - HIDDEN
        for t in linspace(0, math.pi / 2, 10)[1:]:
            cover.append((edge + side * math.sin(t) / abs(sa), top + lift * math.cos(t)))
        cover.append((edge + side / abs(sa), top - 1.2))

    # The free gingiva steepens just enough to meet that surface; the steepness
    # varies continuously with azimuth, so neighbouring sections stay smooth.
    def slope(steep):
        return [
            (p3[0] + k * 0.05, p3[1] - steep * (0.95 * k * 0.05 + 0.75 * (k * 0.05) ** 2))
            for k in range(140)
        ]

    def meet(cone):
        for i in range(len(cone) - 1):
            for j in range(len(cover) - 1):
                t = crossing(cone[i], cone[i + 1], cover[j], cover[j + 1])
                if t is not None:
                    return i, j, t
        return None

    lo, hi = 0.05, 8.0
    for _ in range(14):
        mid = (lo + hi) / 2
        if meet(slope(mid)):
            hi = mid
        else:
            lo = mid
    cone = slope(max(1.0, 1.15 * hi))
    hit = meet(cone)
    assert hit, ("Free gingival slope does not meet the cover", health, a)
    i, j, t = hit
    x = tuple(cone[i][q] + t * (cone[i + 1][q] - cone[i][q]) for q in range(2))
    before = [p for p in cone[: i + 1] if math.dist(p, x) > 0.30]
    after = [p for p in cover[j + 1 :] if math.dist(p, x) > 0.30]
    if not after:
        after = [cover[-1]]
    fillet = bezier(before[-1], x, x, after[0], 12)
    outer = before[:-1] + fillet + after[1:]
    split = len(before) - 1 + len(fillet)
    end = outer[-1]
    y_bot = min(min(bone_top(r) for r in linspace(0.3, end[0], 24)) - 0.7, end[1] - 1.8)
    # Round the turn from the hidden surface down into hard tissue.
    outer += bezier(
        end,
        (end[0] + 0.1, end[1] - 0.3),
        (end[0] - 0.6, end[1] - 0.8),
        (end[0] - 0.6, end[1] - 1.2),
        10,
    )[1:]

    def section(inset):
        deep = 0.15 if inset else 0.0
        lower = 4.0 * inset
        wall_bottom = attach + 0.12 - 2.4 * inset
        wall = [(rt(y) + gap(y) + inset, y) for y in linspace(wall_bottom, wall_top - lower, 40)]
        # Only the visible free gingiva needs a true offset; the hidden remainder
        # simply runs deeper, which avoids folding at its turn into hard tissue.
        chain = (
            offset_chain(outer[:split], inset) + [(r, y - inset - 0.15) for r, y in outer[split:]]
            if inset
            else outer
        )
        r0, y0 = wall[-1]
        if inset:
            # The core begins where the free gingiva is thick enough to hold it,
            # leaving a thicker epithelial cap over the margin.
            k = 0
            while chain[k][1] > y0 + 0.15 or chain[k][0] < r0 + 0.10:
                k += 1
            chain = chain[k:]
            if chain[0][1] < y0 + 0.05:
                # Thin free gingiva: end the core wall below where the core begins.
                wall = [p for p in wall if p[1] < chain[0][1] - 0.2] or wall[:4]
                r0, y0 = wall[-1]
            margin_arc = bezier(
                (r0, y0), (r0, y0 + 0.15), (chain[0][0] - 0.12, chain[0][1] + 0.10), chain[0], 20
            )
        else:
            margin_arc = bezier(
                (r0, y0),
                (r0 - 0.01, margin + 0.06),
                (r0 + 0.20, margin + 0.14),
                chain[0],
                20,
            )
        re, ye = chain[-1]
        embedded = [(re - deep, y) for y in linspace(ye, y_bot - deep, 10)]
        floor = [(r, y_bot - deep) for r in linspace(re - deep, rt(y_bot) - 0.45 - deep, 14)]
        root_y = wall_bottom - 0.22
        along = [(rt(y) - 0.45 - deep, y) for y in linspace(y_bot - deep, root_y, 24)]
        rr = rt(root_y)
        fundus = bezier(
            along[-1],
            (rr + 0.02, root_y),
            (wall[0][0], root_y + 0.05),
            wall[0],
            14,
        )
        parts = [
            (wall, 22),
            (margin_arc, 12),
            (chain, 84),
            (embedded, 6),
            (floor, 10),
            (along, 14),
            (fundus, 8),
        ]
        profile = []
        for points, n in parts:
            profile.extend(resample(points, n)[:-1])
        assert simple_polygon(profile), ("Self-intersecting soft-tissue section", health, a, inset)
        return profile

    # Sulcus or pocket space between the tooth and the gum wall, used to keep the
    # continuous gingiva out of it. It stops just short of the collar wall.
    ys = linspace(attach + 0.05, margin + 0.25, 20)
    sulcus = [(rt(y) - 0.3, y) for y in ys] + [(rt(y) + gap(y) - 0.03, y) for y in ys[::-1]]
    return section(0.0), section(EPITHELIUM), sulcus


def soft_volumes(health, tooth_objects, na=192):
    tooth = surface_tree(tooth_objects)
    rings = {0: [], 1: [], 2: []}
    for i in range(na):
        a = i / na * math.tau
        for k, profile in enumerate(soft_sections(health, a, tooth)):
            rings[k].append(profile)
    out = []
    for k, (name, tid) in enumerate(
        [
            ("Selected-tooth gingival collar", "gingiva"),
            ("Selected-tooth connective core", "connective"),
            ("Selected-tooth sulcus space", "lumen"),
        ]
    ):
        nr = len(rings[k][0])
        assert all(len(p) == nr for p in rings[k])
        pts = []
        for i, profile in enumerate(rings[k]):
            a = i / na * math.tau
            pts.extend((r * math.cos(a), y, r * math.sin(a)) for r, y in profile)
        faces = []
        for i in range(na):
            for j in range(nr):
                faces.append(
                    (
                        i * nr + j,
                        i * nr + (j + 1) % nr,
                        ((i + 1) % na) * nr + (j + 1) % nr,
                        ((i + 1) % na) * nr + j,
                    )
                )
        out.append(make(name, tid, pts, faces))
    return out


def plaque_films(health, tooth_objects):
    """Thin biofilm compartments following the natural crown and root surface.

    Supragingival: cervical crown above the gingival margin, all the way round.
    Subgingival: inside the sulcus on the buccal (cheek) and lingual (tongue)
    sides, stopping short of the selected distal site, whose existing ribbon
    keeps that compartment. These are compartments, not species positions.
    """
    tooth = surface_tree(tooth_objects)
    out = []
    skip = 0.34
    for name, tid, sub in [
        ("Supragingival plaque film", "supragingival", False),
        ("Circumferential subgingival plaque film", "plaque", True),
    ]:
        na = 180 if sub else 192
        angles = (
            [skip + (math.tau - 2 * skip) * i / (na - 1) for i in range(na)]
            if sub
            else [i / na * math.tau for i in range(na)]
        )
        profiles = []
        for a in angles:
            _, margin, attach = soft_landmarks(health, a)
            lo, hi = (attach + 0.15, margin - 0.02) if sub else (margin + 0.02, margin + 1.05)
            ys = linspace(lo, hi, 24)
            last = 4.5
            inner = []
            for y in ys:
                r = outer_radius(tooth, a, y)
                last = r if r is not None and r > 1.2 else last
                inner.append((last + 0.008, y))
            outer = [
                (r + 0.012 + 0.045 * math.sin(math.pi * k / 23) ** 0.5, y)
                for k, (r, y) in enumerate(inner)
            ]
            profiles.append(inner + outer[::-1])
        n = len(profiles[0])
        pts = [
            (r * math.cos(a), y, r * math.sin(a))
            for a, prof in zip(angles, profiles, strict=True)
            for r, y in prof
        ]
        faces = []
        steps = na - 1 if sub else na
        for i in range(steps):
            for j in range(n):
                i2 = (i + 1) % na
                faces.append((i * n + j, i * n + (j + 1) % n, i2 * n + (j + 1) % n, i2 * n + j))
        if sub:
            half = n // 2
            for base in (0, (na - 1) * n):
                for k in range(half - 1):
                    faces.append((base + k, base + k + 1, base + n - 2 - k, base + n - 1 - k))
        out.append(make(name, tid, pts, faces))
    return out


def place(objects, scale, tx):
    """Scale about the 36 origin, then move along the arch to x = tx."""
    tz = bend(tx)
    for o in objects:
        for v in o.data.vertices:
            x, y, z = v.co.x, v.co.z, -v.co.y
            v.co = Vector(coords((x * scale + tx, y * scale, z * scale + tz)))
        o.data.update()


def trim_coronal(roots):
    # Remove the coronal part of a root envelope; the crown covers it.
    bpy.ops.mesh.primitive_cube_add(size=1, location=coords((0, 15.06, 0)))
    cutter = bpy.context.object
    cutter.scale = (40, 40, 30)
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    boolean(roots, cutter)
    discard(cutter)


def above(obj, y):
    # Keep the part of obj above height y.
    bpy.ops.mesh.primitive_cube_add(size=1, location=coords((0, y + 20, 0)))
    box = bpy.context.object
    box.scale = (60, 60, 40)
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    boolean(obj, box, "INTERSECT")
    discard(box)


def support(roots, label):
    cement_outer = combine(
        [expanded(roots, 0.105, label + " cementum envelope", "cementum")],
        label + " cementum envelope",
        "cementum",
        0.08,
    )
    cement = shell_from(roots, cement_outer, label + " cementum", "cementum")
    pdl_outer = combine(
        [expanded(cement_outer, 0.33, label + " socket envelope", "pdl")],
        label + " socket envelope",
        "pdl",
        0.10,
    )
    lining = combine(
        [expanded(pdl_outer, 0.18, label + " socket lining", "bone")],
        label + " socket lining",
        "bone",
        0.10,
    )
    return cement_outer, cement, pdl_outer, lining


def sulcus_cutter(cement_outer, solid, attach, label):
    # Narrow sulcus around a neighbour: the tooth envelope grown slightly and kept
    # above its epithelial attachment. The gingival margin becomes a thin edge.
    env = combine(
        [
            expanded(cement_outer, 0.10, label + " root sulcus", "gingiva"),
            expanded(solid, 0.10, label + " crown sulcus", "gingiva"),
        ],
        label + " sulcus space",
        "gingiva",
        0.07,
    )
    above(env, attach)
    return env


def premolar():
    """Mandibular second premolar (FDI 35) in its own local frame."""
    inner = crown_mesh(True, PREMOLAR, "FDI 35 coronal dentin")
    enamel = crown_mesh(False, PREMOLAR, "FDI 35 enamel")
    solid = crown_mesh(False, PREMOLAR, "FDI 35 crown envelope")
    boolean(enamel, inner)
    rows = []
    for t in linspace(0, 1, 76):
        taper = root_taper(t)
        rows.append(
            (
                0.25 * t * t,
                -2 - 13.2 * t,
                0.05 * math.sin(math.pi * t),
                1.42 * taper + 0.03,
                2.05 * taper + 0.03,
            )
        )

    def concavity(a, t):
        along = smoothstep(t / 0.12) * (1 - smoothstep((t - 0.6) / 0.35))
        return 0.12 * math.exp(-(math.sin(a) ** 2) / 0.22) * along

    roots = combine(
        [
            sphere("FDI 35 root trunk", "dentin", (0, -1.0, 0), (3.25, 2.6, 3.6)),
            rings_volume("FDI 35 root", "dentin", rows, 72, concavity),
        ],
        "FDI 35 root envelope",
        "dentin",
        0.11,
    )
    trim_coronal(roots)
    dentin = combine([inner, clone(roots, "FDI 35 dentin root")], "FDI 35 dentin", "dentin", 0.11)
    parts = [superellipsoid("FDI 35 pulp chamber", "pulp", (0, 1.0, 0), (1.15, 1.0, 1.55), 0.6)]
    for cz, tip in ((0.95, 3.0), (-0.85, 2.3)):
        horn = [
            (
                0,
                1.4 + t * (tip - 1.4),
                cz * (0.4 + 0.6 * t),
                0.45 * (1 - t) ** 0.7 + 0.05,
                0.40 * (1 - t) ** 0.7 + 0.045,
            )
            for t in linspace(0, 1, 30)
        ]
        parts.append(rings_volume("FDI 35 pulp horn", "pulp", horn, 32))
    canal = []
    for y in linspace(0.4, -15.2, 84):
        u = (0.4 - y) / 15.6
        t = max(0.0, (-2 - y) / 13.2)
        r = 0.42 * (1 - u) ** 0.9 + 0.04
        canal.append((0.25 * t * t, y, 0.05 * math.sin(math.pi * t), r, r * 1.35))
    parts.append(rings_volume("FDI 35 root canal", "pulp", canal, 32))
    pulp = combine(parts, "FDI 35 pulp", "pulp", 0.065)
    boolean(dentin, pulp)
    cement_outer, cement, pdl_outer, lining = support(roots, "FDI 35")
    discard(roots)
    return {
        "enamel": enamel,
        "dentin": dentin,
        "pulp": pulp,
        "cement": cement,
        "cement_outer": cement_outer,
        "pdl_outer": pdl_outer,
        "lining": lining,
        "solid": solid,
    }


def furcation_filler():
    # Soft tissue fills any furcation space left above the interradicular bone;
    # the tooth and bone boundaries remove the rest of this box.
    bpy.ops.mesh.primitive_cube_add(size=1, location=coords((0, -5.775, 0)))
    o = bpy.context.object
    o.scale = (3.2, 6.0, 4.85)
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    o.select_set(False)
    return o


def tidy(o):
    # Remove detached Boolean fragments far smaller than the tissue body.
    bm = bmesh.new()
    bm.from_mesh(o.data)
    unseen = set(bm.verts)
    groups = []
    while unseen:
        group = {unseen.pop()}
        todo = list(group)
        while todo:
            v = todo.pop()
            for edge in v.link_edges:
                other = edge.other_vert(v)
                if other in unseen:
                    unseen.remove(other)
                    group.add(other)
                    todo.append(other)
        groups.append(group)
    membership = {v: i for i, g in enumerate(groups) for v in g}
    volumes = [0.0] * len(groups)
    for f in bm.faces:
        volumes[membership[f.verts[0]]] += f.calc_center_median().dot(f.normal) * f.calc_area() / 3
    main = max(abs(v) for v in volumes)
    minor = [i for i, v in enumerate(volumes) if abs(v) < 0.002 * main]
    fraction = sum(abs(volumes[i]) for i in minor) / main
    assert fraction < 0.001, ("Substantial detached soft tissue requires review", o.name, fraction)
    if minor:
        bmesh.ops.delete(bm, geom=[v for i in minor for v in groups[i]], context="VERTS")
    print("SOFT_COMPONENTS", o.name, len(groups), "removed", len(minor), fraction, flush=True)
    o["removed_boolean_fragment_components"] = len(minor)
    o["removed_boolean_fragment_volume_fraction"] = fraction
    bm.to_mesh(o.data)
    bm.free()


def side_of_tooth(y, z=0):
    # Analytic authoring reference for the selected distal pocket compartment.
    neck = (
        4.65 * math.sqrt(max(0, 1 - ((y + 1.2) / 3.2) ** 2 - (z / 3.8) ** 2))
        if -4.4 <= y <= 2
        else 0
    )
    t = max(0, min(1, (-y - 2) / 11.60))
    x, _, cz = root_center(1, t)
    rx, rz = root_radii(1, t)
    rt = x + rx * math.sqrt(max(0, 1 - ((z - cz) / rz) ** 2)) if y <= -2 else 0
    return max(neck, rt) + 0.13


def pocket_ribbon(name, tid, health, offset, thickness, extend=0):
    margin = 0.35 if health == "healthy" else -0.15
    floor = -1.65 if health == "healthy" else -6.20
    na = 31
    nr = 67
    pts = []
    faces = []
    # Paired surfaces, rounded narrow entrance and a smoothly wider interior.
    for outer in [False, True]:
        for j in range(nr):
            t = j / (nr - 1)
            y = margin + (floor - extend - margin) * t
            bulge = (
                0.17 + 0.35 * math.sin(math.pi * t) ** 1.2
                if health == "periodontitis"
                else 0.13 + 0.12 * math.sin(math.pi * t)
            )
            for i in range(na):
                z = -1.18 + i / (na - 1) * 2.36
                width = (
                    offset
                    + (bulge if tid == "lumen" and outer else 0)
                    + (thickness if outer else 0)
                )
                if tid == "epithelium":
                    width += bulge
                pts.append((side_of_tooth(y, z) + width, y, z))
    count = na * nr
    for j in range(nr - 1):
        for i in range(na - 1):
            f = (j * na + i, j * na + i + 1, (j + 1) * na + i + 1, (j + 1) * na + i)
            faces.extend([tuple(reversed(f)), tuple(count + q for q in f)])
    for j in range(nr - 1):
        for i in [0, na - 1]:
            a = j * na + i
            b = (j + 1) * na + i
            faces.append((a, b, count + b, count + a))
    for j in [0, nr - 1]:
        for i in range(na - 1):
            a = j * na + i
            b = a + 1
            faces.append((a, b, count + b, count + a))
    return make(name, tid, pts, faces)


COLORS = {
    "enamel": (0.90, 0.86, 0.76, 1),
    "dentin": (0.78, 0.62, 0.42, 1),
    "pulp": (0.62, 0.22, 0.20, 1),
    "cementum": (0.70, 0.58, 0.40, 1),
    "pdl": (0.55, 0.40, 0.45, 1),
    "gingiva": (0.80, 0.36, 0.32, 1),
    "connective": (0.80, 0.52, 0.46, 1),
    "bone": (0.80, 0.74, 0.62, 1),
    "plaque": (0.55, 0.47, 0.24, 1),
    "supragingival": (0.70, 0.63, 0.40, 1),
    "lumen": (0.16, 0.33, 0.36, 1),
    "epithelium": (0.86, 0.50, 0.44, 1),
    "vessels": (0.58, 0.04, 0.06, 1),
    "nerve": (0.74, 0.52, 0.13, 1),
}
MATS = {}
for tid, color in COLORS.items():
    m = bpy.data.materials.new(tid)
    m.diffuse_color = color
    m.use_nodes = True
    n = m.node_tree.nodes
    bs = n.get("Principled BSDF")
    bs.inputs["Base Color"].default_value = color
    bs.inputs["Roughness"].default_value = (
        0.30 if tid == "enamel" else 0.41 if tid in ["gingiva", "pulp"] else 0.68
    )
    if tid in ["enamel", "gingiva", "epithelium", "pulp"]:
        bs.inputs["Coat Weight"].default_value = 0.24 if tid != "enamel" else 0.48
        bs.inputs["Subsurface Weight"].default_value = 0.045 if tid == "enamel" else 0.09
    noise = n.new("ShaderNodeTexNoise")
    noise.inputs["Scale"].default_value = 16 if tid == "bone" else 60
    noise.inputs["Detail"].default_value = 3
    bump = n.new("ShaderNodeBump")
    bump.inputs["Strength"].default_value = 0.10 if tid == "bone" else 0.035
    bump.inputs["Distance"].default_value = 0.015 if tid == "bone" else 0.004
    m.node_tree.links.new(noise.outputs["Fac"], bump.inputs["Height"])
    m.node_tree.links.new(bump.outputs["Normal"], bs.inputs["Normal"])
    MATS[tid] = m

TRABECULAR = MATS["bone"].copy()
TRABECULAR.name = "Schematic trabecular cut-face material"
nodes = TRABECULAR.node_tree.nodes
links = TRABECULAR.node_tree.links
vor = nodes.new("ShaderNodeTexVoronoi")
vor.feature = "DISTANCE_TO_EDGE"
vor.inputs["Scale"].default_value = 34.0
ramp = nodes.new("ShaderNodeValToRGB")
ramp.color_ramp.elements[0].position = 0.06
ramp.color_ramp.elements[0].color = (0.82, 0.76, 0.62, 1)
ramp.color_ramp.elements[1].position = 0.13
ramp.color_ramp.elements[1].color = (0.42, 0.22, 0.17, 1)
links.new(vor.outputs["Distance"], ramp.inputs[0])
links.new(ramp.outputs["Color"], nodes.get("Principled BSDF").inputs["Base Color"])
# Common tooth geometry is authored independently of the disease presets.
inner_crown = crown_mesh(True)
outer_crown = crown_mesh(False)
crown_solid = crown_mesh(False, MOLAR, "FDI 36 crown envelope")
boolean(outer_crown, inner_crown)
outer_crown.name = "Five-cusp enamel shell"
roots = combine(
    [sphere("Domed root trunk", "dentin", (0, -1.2, 0), (4.65, 3.2, 3.8)), root(-1), root(1)],
    "Root envelope",
    "dentin",
    0.11,
)
# A curved interradicular roof removes the primitive flat trunk underside.
furcation = sphere("Domed furcation authoring cutter", "dentin", (0, -8.7, 0), (1.23, 5.15, 4.80))
boolean(roots, furcation)
discard(furcation)
# Remove the coronal part of the root envelope for the cementum reference.
trim_coronal(roots)
dentin = combine(
    [inner_crown, clone(roots, "Dentin root copy")],
    "Continuous crown and root dentin",
    "dentin",
    0.11,
)
pulp = pulp_body()
boolean(dentin, pulp)
cement_outer = combine(
    [expanded(roots, 0.105, "Cementum external envelope", "cementum")],
    "Cementum external envelope",
    "cementum",
    0.08,
)
cement = shell_from(roots, cement_outer, "Continuous cementum covering", "cementum")
pdl_outer = combine(
    [expanded(cement_outer, 0.33, "Socket and ligament external envelope", "pdl")],
    "Socket and ligament external envelope",
    "pdl",
    0.10,
)
pdl_template = shell_from(cement_outer, pdl_outer, "Periodontal ligament envelope", "pdl")
# A thin cortical lining (lamina dura) surrounds the socket inside the trabecular bone.
socket_lining = combine(
    [expanded(pdl_outer, 0.18, "Socket lining envelope", "bone")],
    "Socket lining envelope",
    "bone",
    0.10,
)

# Neighbouring teeth. 37 reuses the finished 36 parts, scaled; 35 is a premolar.
neighbours = {"35": premolar()}
neighbours["37"] = {
    key: clone(obj, "FDI 37 " + key.replace("_", " "))
    for key, obj in {
        "enamel": outer_crown,
        "dentin": dentin,
        "pulp": pulp,
        "cement": cement,
        "cement_outer": cement_outer,
        "pdl_outer": pdl_outer,
        "lining": socket_lining,
        "solid": crown_solid,
    }.items()
}
for key, (tx, scale, _) in NEIGHBOURS.items():
    n = neighbours[key]
    n["sulcus"] = sulcus_cutter(n["cement_outer"], n["solid"], -1.5, "FDI " + key)
    place(list(n.values()), scale, tx)


def joined(objects, name, tid):
    return combine([clone(o, name, tid) for o in objects], name, tid)


sockets = joined(
    [pdl_outer, *(n["pdl_outer"] for n in neighbours.values())], "All socket envelopes", "pdl"
)
linings = joined(
    [socket_lining, *(n["lining"] for n in neighbours.values())], "All socket linings", "bone"
)

states = {}
for health in ["healthy", "periodontitis"]:
    tooth = [clone(o, o.name, o["tissueId"]) for o in [outer_crown, dentin, pulp, cement]]
    context = []
    for key, n in neighbours.items():
        for part, label in [
            ("enamel", "enamel"),
            ("dentin", "dentin"),
            ("pulp", "pulp"),
            ("cement", "cementum"),
        ]:
            context.append(clone(n[part], "Neighbouring FDI " + key + " " + label))
    bone_whole = bone_tube(health)
    interior_long = bone_tube(health, 0.72)
    boolean(interior_long, linings)
    jaw = clone(bone_whole, "Buccal and lingual cortical plates", "bone")
    boolean(jaw, interior_long)
    boolean(jaw, sockets)
    interior = clone(interior_long, "Trabecular alveolar interior", "bone")
    boolean(interior, bone_whole, "INTERSECT")
    discard(interior_long)
    # No cortical remesh: the plates keep the exact undivided outer surface that
    # the soft tissues are cut against, so gum and bone meet without a seam.
    pdl = clone(pdl_outer, "Periodontal ligament at alveolar socket")
    boolean(pdl, bone_whole, "INTERSECT")
    boolean(pdl, cement_outer)
    for key, n in neighbours.items():
        ligament = clone(n["pdl_outer"], "Neighbouring FDI " + key + " ligament", "pdl")
        boolean(ligament, bone_whole, "INTERSECT")
        boolean(ligament, n["cement_outer"])
        context.append(ligament)
    gum, connective = slab_volumes(health)
    collar, collar_core, sulcus = soft_volumes(health, [outer_crown, cement_outer])
    boolean(gum, collar, "UNION")
    boolean(connective, collar_core, "UNION")
    filler = furcation_filler()
    boolean(connective, filler, "UNION")
    for o in (collar, collar_core, filler):
        discard(o)
    # Partition the intact gingiva first, then carve the same socket, tooth and
    # pocket boundaries from each part. Hidden faces sit inside hard tissue or
    # below the gum surface, so these cuts cross steeply instead of grazing.
    boolean(gum, connective)
    lumen = pocket_ribbon(
        "Fluid lumen with narrow entrance and rounded depth contour", "lumen", health, 0.075, 0.035
    )
    plaque = pocket_ribbon("Tooth-attached plaque compartment", "plaque", health, 0.010, 0.055)
    epi = pocket_ribbon(
        "Pocket lining and junctional attachment", "epithelium", health, 0.075, 0.17, 0.38
    )
    cutters = [bone_whole, cement_outer, crown_solid, sulcus]
    for n in neighbours.values():
        cutters += [n["cement_outer"], n["solid"], n["sulcus"]]
    for part in [gum, connective]:
        for boundary in [*cutters, lumen, epi]:
            boolean(part, boundary)
        tidy(part)
    discard(sulcus)
    films = plaque_films(health, [outer_crown, cement_outer])
    objects = [*tooth, pdl, jaw, interior, gum, connective, plaque, *films, lumen, epi, *context]
    # Small illustrative vascular/neural routes, confined to the authored pulp.
    for tid, r, shift in [("vessels", 0.070, -0.08), ("nerve", 0.045, 0.09)]:
        rows = []
        for j in range(76):
            t = j / 75
            x, y, z = root_center(1, t)
            rows.append((x + shift, y, z, r * (1 - 0.8 * t), r * (1 - 0.8 * t)))
        objects.append(rings_volume("Schematic " + tid + " route", tid, rows, 16))
    discard(bone_whole)
    for o in objects:
        bb = bmesh.new()
        bb.from_mesh(o.data)
        print(
            "PART_CHECK",
            health,
            o.name,
            len(o.data.vertices),
            bb.calc_volume(signed=True),
            sum(not e.is_manifold for e in bb.edges),
            flush=True,
        )
        bb.free()
        o.name = health + " · " + o.name
        o.data.materials.clear()
        o.data.materials.append(TRABECULAR if "Trabecular" in o.name else MATS[o["tissueId"]])
        o.hide_render = health != "periodontitis"
    states[health] = objects
for o in [
    outer_crown,
    crown_solid,
    dentin,
    pulp,
    cement,
    roots,
    cement_outer,
    pdl_outer,
    pdl_template,
    socket_lining,
    sockets,
    linings,
    *(o for n in neighbours.values() for o in n.values()),
]:
    discard(o)


def safe_triangles(m, pos):
    """Exported triangles that keep every edge shared by exactly two faces.

    Blender's n-gon triangulation of Boolean output occasionally picks a
    diagonal that duplicates an existing edge, or the same diagonal as a
    neighbouring n-gon. Such faces are fanned from an added centre point.
    """
    m.calc_loop_triangles()
    by_poly = {}
    for t in m.loop_triangles:
        by_poly.setdefault(t.polygon_index, []).append(tuple(t.vertices))
    used = {tuple(sorted(e.vertices)) for e in m.edges}
    out = []
    for poly in m.polygons:
        tris = by_poly[poly.index]
        vs = list(poly.vertices)
        n = len(vs)
        if n > 3:
            rim = {tuple(sorted((vs[i], vs[(i + 1) % n]))) for i in range(n)}
            diag = {
                tuple(sorted(pair))
                for t in tris
                for pair in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))
            } - rim
            if diag & used:
                centre = len(pos) // 3
                for axis in range(3):
                    pos.append(round(sum(pos[v * 3 + axis] for v in vs) / n, 6))
                tris = [(vs[i], vs[(i + 1) % n], centre) for i in range(n)]
                diag = set()
            used |= diag
        for t in tris:
            out.extend(t)
    return out


export = {
    "format": "marse.pocket-organic-geometry/1",
    "authoring": "Blender " + bpy.app.version_string,
    "units": "illustrative model units",
    "calibration": None,
    "review": "pending",
    "source_sha256": REFERENCE_SCAFFOLD_SHA256,
    "source_role": (
        "Scaffold reference identity only; geometry regenerated from authored "
        "Python rules, not calibrated scans"
    ),
    "cases": {},
    "topology": [],
    "authored_features": {
        "cusps": [{"x": x, "z": z, "height_parameter": a} for x, z, a in CUSPS],
        "orientation": (
            "X mesial (-) to distal (+); Z buccal (+) to lingual (-); Y coronal (+) to apical (-)"
        ),
        "atlas_mesh_reuse": False,
        "pocket_landmarks": {
            "healthy": {"margin": 0.35, "attachment": -1.65},
            "periodontitis": {"margin": -0.15, "attachment": -6.20},
        },
        "layer_widths": "Illustratively expanded for viewing; not physical histometry",
        "segment": {
            "teeth": ["FDI 35", "FDI 36", "FDI 37"],
            "selected": "FDI 36",
            "ends": list(SEGMENT_ENDS),
            "arch_bend": BEND,
            "neighbours": {
                "FDI " + key: {"x": tx, "scale": scale}
                for key, (tx, scale, _) in NEIGHBOURS.items()
            },
            "note": "Authored context segment; neighbours are not reviewed anatomy",
        },
    },
    "authoring_repairs": [],
}
for health, objects in states.items():
    records = []
    for o in objects:
        m = o.data
        pos = []
        for v in m.vertices:
            p = o.matrix_world @ v.co
            pos.extend([round(p.x, 6), round(p.z, 6), round(-p.y, 6)])
        ind = safe_triangles(m, pos)
        assert len(pos) // 3 <= 400000, (o.name, len(pos) // 3)
        # Closed-surface check on the exported triangles themselves.
        uses = {}
        for i in range(0, len(ind), 3):
            for k in range(3):
                edge = tuple(sorted((ind[i + k], ind[i + (k + 1) % 3])))
                uses[edge] = uses.get(edge, 0) + 1
        triangle_edges = sum(c != 2 for c in uses.values())
        records.append(
            {
                "id": o["tissueId"],
                "name": o.get("partName", o.name.split(" · ", 1)[1]),
                "positions": pos,
                "indices": ind,
            }
        )
        bm = bmesh.new()
        bm.from_mesh(m)
        boundary = sum(not e.is_manifold for e in bm.edges) + triangle_edges
        unseen = set(bm.verts)
        components = 0
        while unseen:
            components += 1
            todo = [unseen.pop()]
            while todo:
                v = todo.pop()
                for edge in v.link_edges:
                    other = edge.other_vert(v)
                    if other in unseen:
                        unseen.remove(other)
                        todo.append(other)
        assert len(m.vertices) > 20, (health, o.name, len(m.vertices), boundary)
        assert boundary == 0, (health, o.name, len(m.vertices), boundary)
        repairs = {
            k: o[k]
            for k in o.keys()  # noqa: SIM118 -- Blender ID properties are not iterable
            if k.startswith(("removed_", "filled_", "maximum_", "retained_")) or k == "repair"
        }
        if repairs:
            export["authoring_repairs"].append(
                {"case": health, "part": o.get("partName", o.name.split(" · ", 1)[1]), **repairs}
            )
        bm.free()
        export["topology"].append(
            {
                "case": health,
                "part": o.get("partName", o.name.split(" · ", 1)[1]),
                "nonmanifold_edges": boundary,
                "components": components,
                "triangles": len(ind) // 3,
            }
        )
    export["cases"][health] = records
(OUT / "anatomy.json").write_text(json.dumps(export, separators=(",", ":")) + "\n")
bpy.ops.object.select_all(action="DESELECT")
for o in states["periodontitis"]:
    o.select_set(True)
bpy.ops.export_scene.gltf(
    filepath=str(OUT / "periodontitis.glb"),
    export_format="GLB",
    use_selection=True,
    export_extras=True,
)
scene = bpy.context.scene
world = bpy.data.worlds.new("Dental illustration studio")
scene.world = world
world.use_nodes = True
world.node_tree.nodes["Background"].inputs[0].default_value = (0.017, 0.031, 0.036, 1)
world.node_tree.nodes["Background"].inputs[1].default_value = 0.35


def look(obj, target):
    obj.rotation_euler = (Vector(coords(target)) - obj.location).to_track_quat("-Z", "Y").to_euler()


def light(name, p, energy, color, size):
    data = bpy.data.lights.new(name, "AREA")
    data.energy = energy
    data.color = color
    data.shape = "DISK"
    data.size = size
    o = bpy.data.objects.new(name, data)
    bpy.context.collection.objects.link(o)
    o.location = coords(p)
    look(o, (0, -4, 0))


light("Broad warm key", (-16, 19, 28), 6500, (1, 0.88, 0.76), 20)
light("Soft cool fill", (22, 0, 18), 2500, (0.73, 0.85, 1), 18)
light("Contour light", (-12, 3, -20), 3000, (1, 0.73, 0.55), 16)
bpy.ops.object.camera_add(location=coords((25, 11, 42)))
cam = bpy.context.object
look(cam, (0, -4, 0))
cam.data.type = "ORTHO"
cam.data.ortho_scale = 29
scene.camera = cam
scene.render.engine = "CYCLES"
scene.cycles.samples = 40
scene.cycles.use_denoising = True
scene.render.resolution_x = 1100
scene.render.resolution_y = 1100
scene.render.resolution_percentage = 100
scene.view_settings.view_transform = "AgX"
bpy.context.preferences.filepaths.save_version = 0
bpy.ops.wm.save_as_mainfile(filepath=str(OUT / "pocket-authoring.blend"))
# Full midline reference; browser retains the requested quarter section.
bpy.ops.mesh.primitive_cube_add(size=1, location=coords((0, 0, 25)))
cut = bpy.context.object
cut.scale = (50, 50, 70)
bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
cut.hide_render = True
for rank, o in enumerate(states["periodontitis"]):
    boolean(o, cut)
    for v in o.data.vertices:
        world = o.matrix_world @ v.co
        if abs(world.y) < 0.000001:
            v.co.y -= rank * 0.00005

scene.render.filepath = str(OUT / "blender-cutaway.png")
bpy.ops.render.render(write_still=True)
print("POCKET_ORGANIC_READY", sum(len(o.data.polygons) for o in states["periodontitis"]))
print("NONMANIFOLD_EDGES", sum(x["nonmanifold_edges"] for x in export["topology"]))
