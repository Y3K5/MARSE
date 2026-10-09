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


def occlusal(x, z, grooves=1.0):
    h = 4.30
    for cx, cz, amp in CUSPS:
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
    for mx in (-4.05, 3.95):
        h += 0.20 * math.exp(-(((x - mx) / 0.50) ** 2)) * math.exp(-((z / 2.6) ** 4))
    # Authored central and developmental fissures, not caries or scan data.
    central = z + 0.12 - 0.18 * math.sin(x * 0.9)
    h -= grooves * 0.50 * math.exp(-central * central / 0.045) * math.exp(-x * x / 17)
    for q, sign in GROOVES:
        line = x - q - 0.15 * z
        reach = smoothstep((sign * z + 0.1) / 0.6)
        h -= grooves * 0.30 * math.exp(-line * line / 0.042) * reach
    # Mesial and distal triangular fossae at the ends of the central groove.
    for fx, fz in ((-3.25, 0.0), (2.55, 0.05)):
        h -= grooves * 0.22 * math.exp(-((x - fx) ** 2 + (z - fz) ** 2) / 0.28)
    return h


def cej_height(a):
    # Cervical line curves occlusally on the proximal surfaces.
    return -0.12 + 0.12 * math.cos(2 * a)


def crown_outline(a):
    # Squarer mesial outline and a modest distal taper in occlusal view.
    x = 5.10 * sp(math.cos(a), 0.78) * (1 - 0.046 * math.cos(a))
    z = 4.25 * sp(math.sin(a), 0.83)
    return x, z


def axial_scale(a, t):
    s = math.sin(a)
    # Buccal height of contour in the cervical third, lingual in the middle third;
    # the buccal surface converges more strongly toward the occlusal table.
    peak = 0.40 - 0.13 * s
    top = 0.830 - 0.055 * s
    if t < peak:
        return 0.862 + 0.138 * math.sin(0.5 * math.pi * t / peak)
    return 1 - (1 - top) * ((t - peak) / (1 - peak)) ** 1.7


def crown_mesh(inner=False):
    na = 144
    rows = []
    for j in range(49):
        t = j / 48
        row = []
        for i in range(na):
            a = i / na * math.tau
            ox, oz = crown_outline(a)
            k = axial_scale(a, t)
            # Developmental grooves continue onto the buccal and lingual faces.
            dent = 0
            for q, sign in GROOVES:
                line = ox * k - q - 0.15 * oz * k
                face = smoothstep((sign * math.sin(a) - 0.35) / 0.3)
                dent += 0.07 * math.exp(-line * line / 0.07) * face
            k *= 1 - dent * smoothstep((t - 0.42) / 0.40) / max(math.hypot(ox, oz), 1)
            x = ox * k
            z = oz * k
            cej = cej_height(a)
            y = cej + t * (occlusal(x, z) - cej)
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
            yy = occlusal(xx / 0.846, zz / 0.846, 0.3) - 0.665 if inner else occlusal(xx, zz)
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
    pts.append((0, occlusal(0, 0, 0.3) - 0.665 if inner else occlusal(0, 0), 0))
    for i in range(na):
        faces.extend(
            [
                (bot, i, (i + 1) % na),
                ((len(rows) - 1) * na + i, tip, (len(rows) - 1) * na + (i + 1) % na),
            ]
        )
    return make(
        "Coronal dentin" if inner else "Five-cusp enamel envelope",
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


def crest(a, health):
    loss = (0.5 + 0.5 * math.cos(a)) ** 4 if health == "periodontitis" else 0
    return -3.65 + 0.30 * math.cos(2 * a) - 5.15 * loss


def jaw_volume(health, inset=0):
    na = 128
    nr = 43
    pts = []
    faces = []
    for j in range(nr):
        t = j / (nr - 1)
        for i in range(na):
            a = i / na * math.tau
            rounding = 0.91 + 0.09 * math.sin(t * math.pi * 0.8)
            x = (9.4 - inset) * sp(math.cos(a), 0.46) * rounding
            z = (5.5 - inset) * sp(math.sin(a), 0.62) * rounding + 0.028 * x * x
            top = crest(a, health) - inset * 0.8
            y = -17.0 + inset + t * (top + 17.0 - inset)
            pts.append((x, y, z))
    # The alveolar crest keeps each azimuth's preset height inward toward the
    # tooth, so bone loss stays at the authored distal site instead of sloping
    # down to a single central point on every side.
    # Inside the socket outline the surface blends to the interradicular height.
    # The trabecular interior stays below the cortical top so the furcation keeps
    # a cortical cover; healthy interradicular bone reaches the furcation entrance.
    centre = (-4.2 if health == "periodontitis" else -3.3) - inset * 0.8
    rim = pts[(nr - 1) * na :]
    rings = [(0.86, 1), (0.72, 1), (0.58, 1), (0.44, 1), (0.30, 0.55), (0.18, 0.25)]
    for s, w in rings:
        pts.extend((x * s, centre + w * (y - centre), z * s) for x, y, z in rim)
    rows = nr + len(rings)
    bot = len(pts)
    pts.append((0, -17 + inset, 0))
    top = len(pts)
    pts.append((0, centre, 0))
    for j in range(rows - 1):
        for i in range(na):
            faces.append(
                (j * na + i, (j + 1) * na + i, (j + 1) * na + (i + 1) % na, j * na + (i + 1) % na)
            )
    for i in range(na):
        faces.extend(
            [(bot, i, (i + 1) % na), ((rows - 1) * na + i, top, (rows - 1) * na + (i + 1) % na)]
        )
    return make(
        "Trabecular alveolar interior" if inset else "Buccal and lingual cortical plates",
        "bone",
        pts,
        faces,
    )


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


def top_height(tree, a, r, top=14.0):
    # Height of an upward-facing top surface below (r, a); side faces do not count.
    hit = tree.ray_cast(Vector(coords((r * math.cos(a), top, r * math.sin(a)))), Vector((0, 0, -1)))
    return None if hit[0] is None or hit[1].z < 0.25 else top - hit[3]


def rim(tree, a):
    # Outer edge of the bone's top surface along one azimuth.
    lo, hi = 0.5, 25.0
    for _ in range(32):
        mid = (lo + hi) / 2
        if top_height(tree, a, mid) is None:
            hi = mid
        else:
            lo = mid
    return lo, top_height(tree, a, lo)


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


def soft_landmarks(health, a, rim_y):
    d = (0.5 + 0.5 * math.cos(a)) ** 4 if health == "periodontitis" else 0
    # Scalloped margin: highest at the proximal surfaces, lower mid-buccally and
    # mid-lingually. The selected distal site keeps the preset landmarks exactly.
    margin = 0.35 - 0.30 * math.sin(a) ** 2 - 0.50 * d
    attach = -1.65 - 4.55 * d
    # Rolled basal edge: a smooth line on the cortical plate, kept below the crest.
    base = smin(-7.3 + 0.35 * math.cos(2 * a), rim_y - 2.0, 1.0)
    return d, margin, attach, base


def soft_sections(health, a, tooth, bone):
    """Closed (radius, height) sections of the gingival collar and its connective core.

    Free surfaces are authored against ray-cast tooth and bone surfaces. Parts
    that will later be removed by the socket and cementum boundaries sit well
    inside hard tissue, so every Boolean crossing is steep rather than grazing.
    """
    rim_r, rim_y = rim(bone, a)
    d, margin, attach, base = soft_landmarks(health, a, rim_y)
    tops = [top_height(bone, a, r) for r in linspace(0.3, rim_r - 0.05, 18)]
    y_bot = min(t for t in tops if t is not None) - 0.7
    # Tooth surface radius on a fine height grid. Where a ray passes through the
    # furcation entrance, carry the trunk value downward; the furcation filler
    # added to the connective core covers that region.
    grid = linspace(margin + 0.4, y_bot - 0.3, 220)
    tooth_r = []
    for y in grid:
        r = outer_radius(tooth, a, y)
        tooth_r.append(r if r is not None and r > 1.2 else tooth_r[-1] if tooth_r else 4.6)

    def rt(y):
        f = (grid[0] - y) / (grid[0] - grid[-1]) * (len(grid) - 1)
        i = max(0, min(len(grid) - 2, int(f)))
        f -= i
        return tooth_r[i] + (tooth_r[i + 1] - tooth_r[i]) * f

    def rb(y):
        r = outer_radius(bone, a, min(y, rim_y - 0.02))
        return rim_r if r is None else r

    def gap(y):
        u = max(0.0, min(1.0, (y - attach) / (margin - attach)))
        return 0.11 + (0.03 + 0.14 * d) * math.sin(math.pi * u) ** 0.8

    # Outer surface: a convex free-gingival slope from the margin meets the bone
    # contour offset by the soft-tissue cover, with a fillet at the junction.
    wall_top = margin - 0.18
    rw = rt(wall_top) + gap(wall_top)
    p3 = (rw + 0.48, margin - 0.16)
    # Cover thins slightly toward the outer edge of the crest.
    cover = [
        (r, (top_height(bone, a, r) or rim_y) + COVER_TOP - 0.24 * (r - p3[0]) / (rim_r - p3[0]))
        for r in linspace(p3[0], rim_r, 40)
    ]
    side_y0 = rim_y
    for k in range(1, 13):
        t = k / 12 * math.pi / 2
        cover.append((rim_r + COVER_SIDE * math.sin(t), rim_y + (COVER_TOP - 0.24) * math.cos(t)))
    for y in linspace(side_y0 - 0.1, base + 0.3, 36):
        h = COVER_SIDE - (COVER_SIDE - COVER_BASE) * smoothstep((side_y0 - y) / (side_y0 - base))
        cover.append((rb(y) + h, y))

    # Where the crest is narrow (mid-buccal, mid-lingual) the free gingiva
    # steepens just enough to meet the bone cover. The steepness is a continuous
    # function of azimuth, so neighbouring sections stay smooth.
    def slope(steep):
        return [
            (p3[0] + k * 0.06, p3[1] - steep * (0.95 * k * 0.06 + 0.75 * (k * 0.06) ** 2))
            for k in range(120)
        ]

    def meet(cone):
        for i in range(len(cone) - 1):
            for j in range(len(cover) - 1):
                t = crossing(cone[i], cone[i + 1], cover[j], cover[j + 1])
                if t is not None:
                    return i, j, t
        return None

    lo, hi = 0.2, 8.0
    for _ in range(16):
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
    before = [p for p in cone[: i + 1] if math.dist(p, x) > 0.35]
    after = [p for p in cover[j + 1 :] if math.dist(p, x) > 0.35]
    fillet = bezier(before[-1], x, x, after[0], 12)
    outer = before[:-1] + fillet + after[1:]
    # Rolled basal edge entering the cortical plate steeply.
    end = outer[-1]
    rbb = rb(base)
    outer += bezier(
        end,
        (end[0], base - 0.15),
        (rbb + 0.25, base - 0.45),
        (rbb - 0.30, base - 0.50),
        14,
    )[1:]
    outer.append((rb(base - 0.55) - 0.6, base - 0.55))

    def section(inset):
        deep = 0.15 if inset else 0.0
        lower = 4.0 * inset
        wall_bottom = attach + 0.12 - 2.4 * inset
        wall = [(rt(y) + gap(y) + inset, y) for y in linspace(wall_bottom, wall_top - lower, 40)]
        chain = offset_chain(outer, inset) if inset else outer
        r0, y0 = wall[-1]
        if inset:
            # The core begins where the free gingiva is thick enough to hold it,
            # leaving a thicker epithelial cap over the margin.
            k = 0
            while chain[k][1] > y0 + 0.15 or chain[k][0] < r0 + 0.25:
                k += 1
            chain = chain[k:]
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
        ye = chain[-1][1]
        embedded = [(rb(y) - 0.6 - deep, y) for y in linspace(ye, y_bot - deep, 10)]
        floor = [
            (r, y_bot - deep) for r in linspace(rb(y_bot) - 0.6 - deep, rt(y_bot) - 0.45 - deep, 14)
        ]
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

    return section(0.0), section(EPITHELIUM)


def soft_volumes(health, tooth_objects, bone_object, na=192):
    tooth = surface_tree(tooth_objects)
    bone = surface_tree([bone_object])
    rings = {0: [], 1: []}
    for i in range(na):
        a = i / na * math.tau
        for k, profile in enumerate(soft_sections(health, a, tooth, bone)):
            rings[k].append(profile)
    out = []
    for k, (name, tid) in enumerate(
        [
            ("Continuous scalloped gingival collar", "gingiva"),
            ("Subepithelial connective tissue", "connective"),
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
            _, margin, attach, _ = soft_landmarks(health, a, 0.0)
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
bpy.ops.mesh.primitive_cube_add(size=1, location=coords((0, 15.06, 0)))
cutter = bpy.context.object
cutter.scale = (40, 40, 30)
bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
boolean(roots, cutter)
discard(cutter)
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

states = {}
for health in ["healthy", "periodontitis"]:
    tooth = [clone(o, o.name, o["tissueId"]) for o in [outer_crown, dentin, pulp, cement]]
    jaw = jaw_volume(health)
    interior = jaw_volume(health, 0.72)

    bone_whole = clone(jaw, "Undivided socket boundary")
    boolean(interior, socket_lining)
    boolean(jaw, interior)
    boolean(jaw, pdl_outer)
    # No cortical remesh: the plates keep the exact undivided outer surface that
    # the soft tissues are cut against, so gum and bone meet without a seam.
    pdl = clone(pdl_outer, "Periodontal ligament at alveolar socket")
    boolean(pdl, bone_whole, "INTERSECT")
    boolean(pdl, cement_outer)
    gum, connective = soft_volumes(health, [outer_crown, cement_outer], bone_whole)
    filler = furcation_filler()
    boolean(connective, filler, "UNION")
    discard(filler)
    # Partition the intact collar first, then carve the same socket and pocket
    # boundaries from each part. Embedded faces sit inside hard tissue, so these
    # cuts cross steeply and leave no grazing slivers to remesh.
    boolean(gum, connective)
    lumen = pocket_ribbon(
        "Fluid lumen with narrow entrance and rounded depth contour", "lumen", health, 0.075, 0.035
    )
    plaque = pocket_ribbon("Tooth-attached plaque compartment", "plaque", health, 0.010, 0.055)
    epi = pocket_ribbon(
        "Pocket lining and junctional attachment", "epithelium", health, 0.075, 0.17, 0.38
    )
    for part in [gum, connective]:
        for boundary in [bone_whole, cement_outer, lumen, epi]:
            boolean(part, boundary)
        tidy(part)
    films = plaque_films(health, [outer_crown, cement_outer])
    objects = [*tooth, pdl, jaw, interior, gum, connective, plaque, *films, lumen, epi]
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
    dentin,
    pulp,
    cement,
    roots,
    cement_outer,
    pdl_outer,
    pdl_template,
    socket_lining,
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
        assert len(pos) // 3 <= 65535, (o.name, len(pos) // 3)
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
