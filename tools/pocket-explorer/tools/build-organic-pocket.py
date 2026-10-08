"""Local Blender authoring pipeline. Generic educational geometry, not measured anatomy.
Run with: blender --background --factory-startup --python tools/build-organic-pocket.py
"""

import json
import math
from pathlib import Path

import bmesh
import bpy
from mathutils import Vector

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


CUSPS = [
    (-2.55, 2.05, 1.48),
    (0.70, 2.25, 1.34),
    (3.35, 1.15, 0.95),
    (-2.35, -2.15, 1.76),
    (1.65, -2.0, 1.58),
]


def occlusal(x, z):
    h = 4.70
    for cx, cz, amp in CUSPS:
        h += amp * math.exp(-((x - cx) ** 2 / 1.75 + (z - cz) ** 2 / 1.65))
    # Authored central and developmental fissures, not caries or scan data.
    central = z - 0.18 * math.sin(x * 0.9)
    h -= 0.47 * math.exp(-central * central / 0.050) * math.exp(-x * x / 19)
    for q, sign in [(-0.85, 1), (2.0, 1), (-0.10, -1)]:
        line = x - q - 0.15 * z
        h -= 0.29 * math.exp(-line * line / 0.048) * math.exp(-(((z - sign * 1.0) / 1.8) ** 4))
    h += 0.10 * math.exp(-(((x / 4.7) ** 2 + (z / 3.8) ** 2 - 1) ** 2) / 0.04)
    return h


def crown_mesh(inner=False):
    na = 128
    rows = []
    for j in range(43):
        t = j / 42
        row = []
        for i in range(na):
            a = i / na * math.tau
            wa = 0.862 + 0.136 * math.sin(math.pi * t * 0.87)
            wz = 0.840 + 0.154 * math.sin(math.pi * t * 0.89)
            # Squarer mesial outline, modest distal taper, lingual inclination.
            x = 5.10 * sp(math.cos(a), 0.78) * wa * (1 - 0.046 * math.cos(a))
            z = 4.15 * sp(math.sin(a), 0.83) * wz
            cej = -0.13 + 0.10 * math.cos(2 * a)
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
            yy = occlusal(xx / (0.846 if inner else 1), zz / (0.846 if inner else 1)) - (
                0.665 if inner else 0
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
    pts.append((0, occlusal(0, 0) - (0.665 if inner else 0), 0))
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


def root_radii(side, t):
    return (
        (1.65 if side < 0 else 1.78) * max(0, 1 - t) ** 0.64 + 0.045,
        (2.13 if side < 0 else 1.77) * max(0, 1 - t) ** 0.69 + 0.035,
    )


def rings_volume(name, tid, rows, na=64):
    pts = []
    faces = []
    for x, y, z, rx, rz in rows:
        for i in range(na):
            a = i / na * math.tau
            pts.append((x + rx * math.cos(a), y, z + rz * math.sin(a)))
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
    return rings_volume("Mesial root" if side < 0 else "Distal root", "dentin", rows)


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


def pulp_body():
    objs = [sphere("Pulp chamber roof and floor", "pulp", (0, 1.90, 0), (2.20, 1.37, 1.60))]
    for cx, cz, amp in CUSPS:
        rows = []
        for j in range(30):
            t = j / 29
            rows.append(
                (
                    cx * 0.63 * t,
                    2.35 + t * (1.7 + amp * 0.20),
                    cz * 0.62 * t,
                    0.63 * (1 - t) ** 0.64 + 0.045,
                    0.50 * (1 - t) ** 0.64 + 0.04,
                )
            )
        objs.append(rings_volume("Pulp horn", "pulp", rows, 32))
    for side, shift in [(-1, -0.55), (-1, 0.55), (1, 0)]:
        rows = []
        for j in range(85):
            u = j / 84
            if u < 0.22:
                t = u / 0.22
                x = side * 0.5 + (side * 2.03 - side * 0.5) * smoothstep(t)
                y = 1.9 - 3.9 * t
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
    bot = len(pts)
    pts.append((0, -17 + inset, 0))
    top = len(pts)
    pts.append((0, -7.25 if health == "periodontitis" else -3.8 - inset, 0))
    for j in range(nr - 1):
        for i in range(na):
            faces.append(
                (j * na + i, (j + 1) * na + i, (j + 1) * na + (i + 1) % na, j * na + (i + 1) % na)
            )
    for i in range(na):
        faces.extend(
            [(bot, i, (i + 1) % na), ((nr - 1) * na + i, top, (nr - 1) * na + (i + 1) % na)]
        )
    return make(
        "Trabecular alveolar interior" if inset else "Buccal and lingual cortical plates",
        "bone",
        pts,
        faces,
    )


def tooth_radius(a, y):
    dx, dz = math.cos(a), math.sin(a)
    best = 0
    q = 1 - ((y + 1.2) / 3.2) ** 2
    if q > 0:
        best = math.sqrt(q / (dx * dx / 4.65**2 + dz * dz / 3.8**2))
    if y <= -2:
        for side in [-1, 1]:
            t = (-y - 2) / (12.45 if side < 0 else 11.6)
            if t < 0 or t > 1:
                continue
            cx, _, cz = root_center(side, t)
            rx, rz = root_radii(side, t)
            aa = dx * dx / rx**2 + dz * dz / rz**2
            bb = -2 * (dx * cx / rx**2 + dz * cz / rz**2)
            cc = cx * cx / rx**2 + cz * cz / rz**2 - 1
            disc = bb * bb - 4 * aa * cc
            if disc >= 0:
                best = max(best, (-bb + math.sqrt(disc)) / (2 * aa))
    return best


def catmull(a, b, c, d, t):
    return 0.5 * (
        (2 * b)
        + (-a + c) * t
        + (2 * a - 5 * b + 4 * c - d) * t * t
        + (-a + 3 * b - 3 * c + d) * t * t * t
    )


def gum_volume(health, core=False):
    na = 128
    nr = 96
    pts = []
    faces = []
    for j in range(nr):
        u = j / nr * 10
        k = int(u)
        t = u - k
        for i in range(na):
            a = i / na * math.tau
            d = (0.5 + 0.5 * math.cos(a)) ** 4 if health == "periodontitis" else 0
            margin = 0.35 - 0.50 * d + 0.12 * math.cos(2 * a)
            attach = -1.65 - 4.55 * d
            base = crest(a, health) - 3.0
            keys = [
                (4.43, 3.52, margin),
                (4.58, 3.68, margin + 0.14),
                (5.12, 3.99, margin - 0.26),
                (6.23, 4.58, margin - 1.25),
                (9.53, 5.62, base + 1.9),
                (9.56, 5.69, base),
                (8.60, 4.78, base - 0.1),
                (6.50, 4.05, base + 0.80),
                (4.7, 3.63, attach - 0.50),
                (4.61, 3.64, attach + 0.22),
            ]
            for q, clearance in [(7, 0.35), (8, 0.19), (9, 0.24)]:
                yy = keys[q][2]
                r = tooth_radius(a, yy) + clearance
                keys[q] = (r * abs(math.cos(a)) ** 0.35, r * abs(math.sin(a)) ** 0.26, yy)
            ix = [(k - 1) % 10, k % 10, (k + 1) % 10, (k + 2) % 10]
            xr = catmull(*(keys[q][0] for q in ix), t)
            zr = catmull(*(keys[q][1] for q in ix), t)
            y = catmull(*(keys[q][2] for q in ix), t)
            if core:
                # Inset toward the local cross-section centre: outer surface
                # shrinks while the tooth-facing hole widens, preserving a skin.
                cx = sum(q[0] for q in keys) / 10
                cz = sum(q[1] for q in keys) / 10
                cy = sum(q[2] for q in keys) / 10
                xr = cx + (xr - cx) * 0.83
                zr = cz + (zr - cz) * 0.83
                y = cy + (y - cy) * 0.86
            x = xr * sp(math.cos(a), 0.65)
            z = zr * sp(math.sin(a), 0.74) + 0.018 * x * x
            pts.append((x, y, z))
    for j in range(nr):
        for i in range(na):
            faces.append(
                (
                    j * na + i,
                    j * na + (i + 1) % na,
                    ((j + 1) % nr) * na + (i + 1) % na,
                    ((j + 1) % nr) * na + i,
                )
            )
    return make(
        "Subepithelial connective tissue" if core else "Continuous scalloped gingival collar",
        "connective" if core else "gingiva",
        pts,
        faces,
    )


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
    "enamel": (0.93, 0.88, 0.76, 1),
    "dentin": (0.69, 0.49, 0.25, 1),
    "pulp": (0.38, 0.07, 0.08, 1),
    "cementum": (0.65, 0.55, 0.37, 1),
    "pdl": (0.49, 0.42, 0.57, 1),
    "gingiva": (0.47, 0.16, 0.17, 1),
    "connective": (0.48, 0.23, 0.23, 1),
    "bone": (0.48, 0.39, 0.28, 1),
    "plaque": (0.43, 0.37, 0.16, 1),
    "lumen": (0.13, 0.33, 0.34, 1),
    "epithelium": (0.77, 0.37, 0.30, 1),
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
vor.inputs["Scale"].default_value = 24.0
ramp = nodes.new("ShaderNodeValToRGB")
ramp.color_ramp.elements[0].position = 0.035
ramp.color_ramp.elements[0].color = (0.23, 0.16, 0.09, 1)
ramp.color_ramp.elements[1].position = 0.10
ramp.color_ramp.elements[1].color = (0.66, 0.53, 0.35, 1)
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

states = {}
for health in ["healthy", "periodontitis"]:
    tooth = [clone(o, o.name, o["tissueId"]) for o in [outer_crown, dentin, pulp, cement]]
    jaw = jaw_volume(health)
    interior = jaw_volume(health, 0.72)

    bone_whole = clone(jaw, "Undivided socket boundary")
    boolean(jaw, interior)
    boolean(jaw, pdl_outer)
    boolean(interior, pdl_outer)
    jaw = combine([jaw], "Buccal and lingual cortical plates", "bone", 0.17)
    jaw["repair"] = "Voxel resolution of tiny socket/plate Boolean sliver faces"
    pdl = clone(pdl_outer, "Periodontal ligament at alveolar socket")
    boolean(pdl, bone_whole, "INTERSECT")
    boolean(pdl, cement_outer)
    gum = combine([gum_volume(health)], "Continuous scalloped gingival collar", "gingiva", 0.12)
    connective = combine(
        [gum_volume(health, True)], "Subepithelial connective tissue", "connective", 0.12
    )
    # Partition the intact collar first, then carve the same socket and pocket
    # boundaries from each part. Avoid re-subtracting independently clipped,
    # coincident compartment faces.
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
    # Re-resolve thin CSG slivers as authoring repairs, then retain the common
    # boundary cuts. This is geometric cleanup, never histometric evidence.
    for part in [gum, connective]:
        bm = bmesh.new()
        bm.from_mesh(part.data)
        bad = sum(not e.is_manifold for e in bm.edges)
        bm.free()
        # Regularize thin surfaces even when edge topology is already closed.
        if part in [gum, connective]:
            clean = combine([part], part.name, part["tissueId"], 0.10 if part == gum else 0.12)
            clean["repair"] = "Voxel repair of " + str(bad) + " Boolean nonmanifold edges"
    # Final contact cut uses the rendered cortical/core surfaces, avoiding
    # overlap introduced by the cortical and soft-tissue remeshing passes.
    for soft in [gum, connective]:
        boolean(soft, jaw)
        boolean(soft, interior)
        soft["repair"] = soft.get("repair", "") + "; final contact cut against rendered bone"
    objects = [*tooth, pdl, jaw, interior, gum, connective, plaque, lumen, epi]
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
for o in [outer_crown, dentin, pulp, cement, roots, cement_outer, pdl_outer, pdl_template]:
    discard(o)

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
        m.calc_loop_triangles()
        pos = []
        ind = []
        for v in m.vertices:
            p = o.matrix_world @ v.co
            pos.extend([round(p.x, 6), round(p.z, 6), round(-p.y, 6)])
        for t in m.loop_triangles:
            ind.extend(t.vertices)
        assert len(m.vertices) <= 65535, (o.name, len(m.vertices))
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
        boundary = sum(not e.is_manifold for e in bm.edges)
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
