"""Pack a selected public BodyParts3D subset without altering source coordinates.

Usage: python3 tools/build-atlas-context.py SOURCE_DIRECTORY
The source directory holds the official PART-OF and IS-A 4.0 archives. No
network requests, source execution, or patient data are used by this builder.
"""

import base64
import hashlib
import json
import re
import struct
import sys
import zipfile
from pathlib import Path

source = Path(sys.argv[1])
root = Path(__file__).resolve().parents[1]
archives = [source / f"{tree}_BP3D_4.0_obj_99.zip" for tree in ("partof", "isa")]
chosen = {}
palatal_elements = {}
for line in (source / "partof_element_parts.txt").read_text().splitlines():
    columns = line.split("\t")
    if len(columns) == 3 and columns[0] in {"FMA55021", "FMA55022"}:
        palatal_elements[columns[2]] = dict(concept=columns[0], name=columns[1])
bones = {
    "mandible",
    "right maxilla",
    "left maxilla",
    "right palatine bone",
    "left palatine bone",
    "frontal bone",
    "sphenoid bone",
    "right temporal bone",
    "left temporal bone",
    "right zygomatic bone",
    "left zygomatic bone",
    "right nasal bone",
    "left nasal bone",
    "ethmoid",
    "occipital bone",
    "right parietal bone",
    "left parietal bone",
}
soft = {
    "tongue": "tongue",
    "gingiva of lower jaw": "gingiva",
    "gingiva of upper jaw": "gingiva",
    "uvula": "palate",
    "soft palate": "palate",
    "left submandibular gland": "connective",
    "right submandibular gland": "connective",
    "left sublingual gland": "connective",
    "right sublingual gland": "connective",
}
vessels = {
    f"{side} {name}": style
    for side in ("left", "right")
    for name, style in (
        ("common carotid artery", "artery"),
        ("internal carotid artery", "artery"),
        ("internal jugular vein", "vein"),
    )
}
for archive in archives:
    with zipfile.ZipFile(archive) as z:
        for member in z.namelist():
            if not member.endswith(".obj"):
                continue
            data = z.read(member)
            text = data.decode("utf-8")
            match = re.search(r"^# English name\s*:\s*(.+)$", text, re.M)
            if not match:
                continue
            name = match[1].strip().lower()
            fid = Path(member).stem
            tooth = bool(
                re.fullmatch(
                    r"(left|right) (upper|lower) (?:.* )?secondary "
                    r"(molar|premolar|incisor|canine) tooth",
                    name,
                )
            )
            if (
                name not in bones | set(soft) | set(vessels)
                and not tooth
                and fid not in palatal_elements
            ):
                continue
            if fid in chosen:
                # PART-OF is the selected identity; IS-A may have different
                # representation/header bytes. Record that alternative separately.
                chosen[fid]["alternate_source"] = dict(
                    archive=archive.name, member=member, sha256=hashlib.sha256(data).hexdigest()
                )
                continue
            positions, indices = [], []
            for line in text.splitlines():
                if line.startswith("v "):
                    xyz = [float(x) for x in line.split()[1:4]]
                    if len(xyz) != 3 or any(abs(x) > 1e5 for x in xyz):
                        raise ValueError("Invalid atlas vertex")
                    positions.extend(xyz)
                elif line.startswith("f "):
                    face = [int(x.split("/")[0]) - 1 for x in line.split()[1:]]
                    for j in range(1, len(face) - 1):
                        indices.extend((face[0], face[j], face[j + 1]))
            if (
                not positions
                or not indices
                or min(indices) < 0
                or max(indices) * 3 >= len(positions)
            ):
                raise ValueError("Invalid source indices")
            concept = re.search(r"^# Concept ID\s*:\s*(.+)$", text, re.M)[1].strip()
            jaw = "lower" if name == "mandible" or "lower" in name or name == "tongue" else "upper"
            fdi = None
            if tooth:
                quadrant = (
                    (2 if name.startswith("left") else 1)
                    if "upper" in name
                    else (3 if name.startswith("left") else 4)
                )
                number = (
                    3
                    if "canine" in name
                    else (1 if "central" in name else 2)
                    if "incisor" in name
                    else (4 if "first" in name else 5)
                    if "premolar" in name
                    else (6 if "first" in name else 7)
                )
                fdi = quadrant * 10 + number
            tissue = (
                "enamel"
                if tooth
                else "bone"
                if name in bones
                else "vessels"
                if name in vessels
                else "palate"
                if fid in palatal_elements
                else soft[name]
            )
            chosen[fid] = dict(
                element=fid,
                concept=concept,
                name=name,
                id=tissue,
                fdi=fdi,
                jaw=jaw,
                appearance=vessels.get(name),
                compound=palatal_elements.get(fid),
                scope="mouth"
                if tooth
                or fid in palatal_elements
                or name in soft
                or name in vessels
                or name
                in {
                    "mandible",
                    "right maxilla",
                    "left maxilla",
                    "right palatine bone",
                    "left palatine bone",
                }
                else "face",
                positions=base64.b64encode(
                    struct.pack("<" + "f" * len(positions), *positions)
                ).decode(),
                indices=base64.b64encode(struct.pack("<" + "I" * len(indices), *indices)).decode(),
                triangles=len(indices) // 3,
                source_archive=archive.name,
                source_member=member,
                source_sha256=hashlib.sha256(data).hexdigest(),
            )
if {p["fdi"] for p in chosen.values() if p["fdi"]} != {
    q * 10 + n for q in range(1, 5) for n in range(1, 8)
}:
    raise ValueError("Expected 28 uniquely identified permanent teeth")
payload = dict(
    format="marse.bodyparts3d-context/1",
    source="https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html",
    release="4.0 (2013); official archive license updated 2025-02-27",
    license="CC-BY-4.0",
    attribution=(
        "BodyParts3D, © The Database Center for Life Science licensed under "
        "CC Attribution 4.0 International"
    ),
    license_source="https://dbarchive.biosciencedbc.jp/en/bodyparts3d/lic.html",
    historical_obj_header_license=(
        "CC BY-SA 2.1 Japan; preserved source headers, superseded on official archive license page"
    ),
    specimen="Adult male reference atlas; not an average anatomy or patient reconstruction",
    units="mm in original OBJ Bounds(mm) headers",
    axes="Source +x subject left, +z superior, -y anterior; viewer x,z-1460,-y-140",
    precision=(
        "Original decimal OBJ positions packed float32; indices uint32 little-endian; "
        "source hashes refer to original OBJ bytes"
    ),
    archives=[
        dict(name=f.name, sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in archives
    ],
    changes=(
        "Subset selection, triangulation of polygon faces, float32 packing; "
        "rigid viewer reorientation. Displayed lower jaw opening is separately authored."
    ),
    parts=sorted(chosen.values(), key=lambda p: p["element"]),
)
(root / "src/pocket-atlas-input.json").write_text(json.dumps(payload, separators=(",", ":")) + "\n")
print(
    json.dumps(
        dict(parts=len(chosen), triangles=sum(p["triangles"] for p in chosen.values()), teeth=28)
    )
)
