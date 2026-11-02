"""Bind retained full-coordinate atom/topology context to an existing local MARSE fold.

Requires the already-installed structural parser (Biopython). Never downloads,
reconstructs atoms, creates a study or executes scientific analyses.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import re
import tempfile
from pathlib import Path

MAX_BYTES = 64 * 1024 * 1024
METALS = frozenset(
    ["NA", "MG", "K", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG", "SR", "CS", "BA"]
)


def digest(raw: bytes) -> str:
    """Return a source-byte pin."""
    return hashlib.sha256(raw).hexdigest()


def checked(path: str, expected: str) -> bytes:
    """Read once, then verify exactly those bytes."""
    raw = Path(path).read_bytes()
    if len(raw) > MAX_BYTES or digest(raw) != expected:
        raise ValueError("Atomic source checksum or size differs")
    return raw


def parse_pdb(raw: bytes) -> tuple[list[dict], dict]:
    """Read first-model primary/A conformer atoms; preserve author identities."""
    atoms, omitted, seen = [], {"alternate": 0, "zero_occupancy": 0}, set()
    for line in raw.decode("utf-8").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        alt = line[16:17].strip()
        if alt not in ("", "A"):
            omitted["alternate"] += 1
            continue
        occupancy = float(line[54:60]) if line[54:60].strip() else None
        if occupancy is not None and occupancy <= 0:
            omitted["zero_occupancy"] += 1
            continue
        atom = {
            "i": len(atoms),
            "serial": int(line[6:11]),
            "name": line[12:16].strip(),
            "element": line[76:78].strip().upper(),
            "chain": line[21:22].strip() or "_",
            "residue": int(line[22:26]),
            "insertion": line[26:27].strip(),
            "resname": line[17:20].strip(),
            "altloc": alt,
            "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            "occupancy": occupancy,
            "b": float(line[60:66]) if line[60:66].strip() else None,
            "group": line[:6].strip(),
        }
        key = (atom["chain"], atom["residue"], atom["insertion"], atom["resname"], atom["name"])
        if (
            key in seen
            or not atom["name"]
            or not atom["resname"]
            or not atom["element"]
            or not all(math.isfinite(x) for x in atom["xyz"])
            or any(atom[k] is not None and not math.isfinite(atom[k]) for k in ("occupancy", "b"))
        ):
            raise ValueError("Ambiguous primary atom identity or incomplete element/coordinates")
        seen.add(key)
        atoms.append(atom)
    if not atoms:
        raise ValueError("No retained primary atom coordinates")
    return atoms, omitted


def cif_rows(doc: dict, category: str) -> list[dict]:
    """Read a parser-qualified mmCIF category without conflating missing tokens."""
    columns = {k.split(".", 1)[1]: v for k, v in doc.items() if k.startswith(category + ".")}
    if not columns:
        return []
    lengths = {len(v) for v in columns.values()}
    if len(lengths) != 1:
        raise ValueError("Inconsistent mmCIF category")
    return [
        dict(zip(columns, values, strict=True)) for values in zip(*columns.values(), strict=True)
    ]


def normalized(value: str | None) -> str:
    """Keep absent CIF fields absent."""
    return "" if value in (None, "?", ".") else str(value)


def residue_key(atom: dict) -> tuple:
    """Identify a residue in the retained author frame."""
    return atom["chain"], atom["residue"], atom["insertion"], atom["resname"]


def qualify_atoms(atoms: list[dict], doc: dict, shape: dict, ledger: dict) -> None:
    """Join every shown atom to source CIF and every backbone anchor to the current fold."""
    rows = cif_rows(doc, "_atom_site")
    first = rows[0].get("pdbx_PDB_model_num", "1") if rows else "1"
    source = {}
    for row in rows:
        if row.get("pdbx_PDB_model_num", "1") != first:
            continue
        key = (
            normalized(row.get("auth_asym_id")) or "_",
            int(row.get("auth_seq_id") or row["label_seq_id"]),
            normalized(row.get("pdbx_PDB_ins_code")),
            normalized(row.get("auth_comp_id") or row.get("label_comp_id")),
            normalized(row.get("auth_atom_id") or row.get("label_atom_id")),
            normalized(row.get("label_alt_id")),
        )
        source.setdefault(key, []).append(row)
    for atom in atoms:
        key = (*residue_key(atom), atom["name"], atom["altloc"])
        matches = source.get(key, [])
        if len(matches) != 1:
            raise ValueError("View/source atom identity is missing or ambiguous")
        row = matches[0]
        xyz = [float(row["Cartn_" + k]) for k in "xyz"]
        if math.dist(xyz, atom["xyz"]) > 0.002 or row["type_symbol"].upper() != atom["element"]:
            raise ValueError("View/source atom coordinates or element differ")
        seq = normalized(row.get("label_seq_id"))
        atom.update(
            label_seq=int(seq) if seq else None,
            label_chain=row["label_asym_id"],
            entity=row.get("label_entity_id"),
            charge=normalized(row.get("pdbx_formal_charge")) or None,
            role=("backbone" if atom["name"] in ("N", "CA", "C", "O", "OXT") else "side chain")
            if atom["group"] == "ATOM"
            else "non-polymer component",
        )
    qualify_backbone(atoms, shape, ledger)


def qualify_backbone(atoms: list[dict], shape: dict, ledger: dict) -> None:
    """Require the exact saved C-alpha/name/chain/insertion coordinate join."""
    observed = {
        (a["chain"], a["residue"], a["insertion"]): a
        for a in atoms
        if a["group"] == "ATOM" and a["name"] == "CA"
    }
    count = 0
    for chain, names in zip(shape["chains"], ledger["chains"], strict=True):
        if chain["chain"] != names["chain"]:
            raise ValueError("Atomic ledger chain differs")
        for point, (number, name, insertion) in zip(chain["ca"], names["residues"], strict=True):
            atom = observed.get((chain["chain"], number, insertion))
            if not atom or atom["resname"] != name or math.dist(atom["xyz"], point[1:]) > 0.00001:
                raise ValueError("Atom/backbone residue/name/coordinate correspondence differs")
            count += 1
    if count != len(observed):
        raise ValueError("Atom/backbone coverage differs")


def qualify_pdb_prediction(
    raw: bytes, atoms: list[dict], shape: dict, ledger: dict, protein: dict, context: dict
) -> dict:
    """Admit a retained complete single-chain prediction, never a homolog or complex."""
    aa = dict(
        zip(
            [
                "ALA",
                "ARG",
                "ASN",
                "ASP",
                "CYS",
                "GLN",
                "GLU",
                "GLY",
                "HIS",
                "ILE",
                "LEU",
                "LYS",
                "MET",
                "PHE",
                "PRO",
                "SER",
                "THR",
                "TRP",
                "TYR",
                "VAL",
            ],
            "ARNDCQEGHILKMFPSTWYV",
            strict=True,
        )
    )
    model = re.fullmatch(r"AF-([A-Z0-9]+)-F([1-9][0-9]*)", shape.get("id", ""))
    if (
        not model
        or shape.get("kind") != "prediction"
        or shape.get("admission") != "prediction_sequence_matched"
        or context.get("kind") != "prediction"
        or context.get("id") != shape["id"]
    ):
        raise ValueError("Retained PDB prediction model identity is not qualified")
    lines = raw.decode().splitlines()
    title = " ".join(line[10:].strip() for line in lines if line.startswith("TITLE "))
    if "ALPHAFOLD" not in title or "PREDICTION" not in title or f"({model[1]})" not in title:
        raise ValueError("Prediction header accession differs from retained model")
    anchors = [a for a in atoms if a["group"] == "ATOM" and a["name"] == "CA"]
    sequence = "".join(aa.get(a["resname"], "?") for a in anchors)
    reference = protein.get("reference", {})
    selected = protein.get("subtractive_record", reference)
    seq_hash = digest(sequence.encode())
    if (
        "?" in sequence
        or sequence != reference.get("sequence")
        or len(sequence) != reference.get("length")
        or seq_hash != reference.get("sequence_sha256")
        or seq_hash != selected.get("sequence_sha256")
    ):
        raise ValueError("Prediction full sequence differs from exact retained reference")
    chains = {a["chain"] for a in atoms}
    if (
        len(chains) != 1
        or any(a["group"] != "ATOM" or a["insertion"] for a in atoms)
        or [a["residue"] for a in anchors] != list(range(1, len(sequence) + 1))
    ):
        raise ValueError(
            "PDB prediction needs complete unambiguous single-chain reference numbering"
        )
    chain = anchors[0]["chain"]
    dbref = [line for line in lines if line.startswith("DBREF ")]
    # Retained AlphaFold PDBs extend the DBREF accession field for ten-character IDs.
    # The adapter already refuses insertions; require the complete literal token identity.
    tokens = dbref[0].split() if len(dbref) == 1 else []
    if (
        len(tokens) != 10
        or tokens[2] != chain
        or tokens[5] != "UNP"
        or tokens[6] != model[1]
        or [int(tokens[i]) for i in (3, 4, 8, 9)] != [1, len(sequence), 1, len(sequence)]
    ):
        raise ValueError("Prediction DBREF identity or numbering differs")
    seqres = [line for line in lines if line.startswith("SEQRES")]
    declared = [name for line in seqres for name in line[19:70].split()]
    if any(
        line[11:12] != chain or int(line[13:17]) != len(sequence) for line in seqres
    ) or declared != [a["resname"] for a in anchors]:
        raise ValueError("Prediction SEQRES differs from resolved reference sequence")
    residues = {a["residue"]: a["resname"] for a in anchors}
    if any(
        not {"N", "CA", "C", "O"}.issubset({a["name"] for a in atoms if a["residue"] == n})
        for n in residues
    ):
        raise ValueError("Prediction lacks retained full backbone atoms; no atoms reconstructed")
    for a in atoms:
        if residues.get(a["residue"]) != a["resname"]:
            raise ValueError("Prediction atom residue differs from reference")
        a.update(
            label_seq=None,
            label_chain=None,
            entity=None,
            charge=None,
            polymer_position=a["residue"],
            polymer_chain=a["chain"],
            role="backbone" if a["name"] in ("N", "CA", "C", "O", "OXT") else "side chain",
        )
    qualify_backbone(atoms, shape, ledger)
    return {
        "sequence_sha256": seq_hash,
        "source_accession": model[1],
        "target_accession": protein.get("accession"),
        "target_strain": selected.get("strain"),
        "source_organism_record": " ".join(
            line[10:].strip() for line in lines if line.startswith("SOURCE")
        ),
        "polymer_numbering": (
            "Full reference sequence, SEQRES, DBREF and author numbering agree; "
            "no mmCIF label identifiers supplied"
        ),
        "limit": (
            "Exact retained sequence/model correspondence; source strain phenotype "
            "and interaction sites are not transferred"
        ),
    }


def build_topology(
    atoms: list[dict], doc: dict, fallback: dict, raw_hash: str, fallback_hash: str
) -> tuple[list[dict], dict]:
    """Use retained component definitions, qualified polymer links and explicit connections."""
    groups, bonds, held, seen = {}, [], [], set()
    for atom in atoms:
        groups.setdefault(residue_key(atom), {})[atom["name"]] = atom
    own_defs, fallback_defs = {}, {}
    for target, source in [(own_defs, doc), (fallback_defs, fallback)]:
        for row in cif_rows(source, "_chem_comp_bond"):
            target.setdefault(row["comp_id"], []).append(row)
    components = {
        r["id"]: {"name": r.get("name"), "type": r.get("type")} for r in cif_rows(doc, "_chem_comp")
    }

    def add(a: dict, b: dict, order: str, kind: str, evidence: str, source_hash: str) -> None:
        key = tuple(sorted((a["i"], b["i"])))
        if key[0] == key[1] or key in seen:
            return
        seen.add(key)
        bonds.append(
            {
                "a": key[0],
                "b": key[1],
                "order": order,
                "kind": kind,
                "evidence": evidence,
                "source_sha256": source_hash,
                "distance": round(math.dist(a["xyz"], b["xyz"]), 5),
            }
        )

    for key, members in groups.items():
        comp = key[3]
        definitions = own_defs.get(comp) or fallback_defs.get(comp) or []
        source_hash = raw_hash if comp in own_defs else fallback_hash
        if not definitions and len(members) > 1:
            held.append({"residue": list(key), "reason": "No retained component topology"})
        for atom in members.values():
            atom["component"] = components.get(comp, {"name": comp, "type": "unknown"})
        for row in definitions:
            a, b = members.get(row["atom_id_1"]), members.get(row["atom_id_2"])
            if not a or not b:
                continue  # Missing atoms, including hydrogens, are never created.
            metal = a["element"] in METALS or b["element"] in METALS
            distance = math.dist(a["xyz"], b["xyz"])
            if not 0.4 < distance < (3.5 if metal else 2.6):
                held.append(
                    {
                        "residue": list(key),
                        "reason": "Component bond geometry held",
                        "atoms": [a["name"], b["name"]],
                    }
                )
                continue
            order = (
                "arom"
                if row.get("pdbx_aromatic_flag", "").upper() == "Y"
                else row["value_order"].lower()
            )
            add(
                a,
                b,
                order,
                "component_coordination" if metal else "component_template",
                "Retained chemical-component topology; not an interaction measurement",
                source_hash,
            )
    sequence = {}
    for members in groups.values():
        atom = next(iter(members.values()))
        position = atom.get("polymer_position", atom["label_seq"])
        if atom["group"] == "ATOM" and position is not None:
            sequence[(atom.get("polymer_chain", atom["label_chain"]), position)] = members
    for (chain, position), members in sequence.items():
        a, b = members.get("C"), sequence.get((chain, position + 1), {}).get("N")
        if a and b and 1.0 <= math.dist(a["xyz"], b["xyz"]) <= 1.8:
            add(
                a,
                b,
                "sing",
                "polymer_link",
                (
                    "Qualified consecutive polymer positions plus source-coordinate geometry; "
                    "conventional peptide topology"
                ),
                raw_hash,
            )
    for row in cif_rows(doc, "_struct_conn"):
        selected = []
        if any(
            normalized(row.get("ptnr" + str(i) + "_symmetry")) not in ("", "1_555") for i in (1, 2)
        ):
            held.append(
                {"connection": row.get("id"), "reason": "Symmetry expansion is not admitted"}
            )
            continue
        for i in (1, 2):
            prefix = "ptnr" + str(i)
            candidates = [
                a
                for a in atoms
                if a["chain"] == (normalized(row.get(prefix + "_auth_asym_id")) or "_")
                and str(a["residue"]) == normalized(row.get(prefix + "_auth_seq_id"))
                and a["name"] == row.get(prefix + "_label_atom_id")
                and a["insertion"] == normalized(row.get("pdbx_" + prefix + "_PDB_ins_code"))
                and (
                    not normalized(row.get("pdbx_" + prefix + "_label_alt_id"))
                    or a["altloc"] == normalized(row.get("pdbx_" + prefix + "_label_alt_id"))
                )
            ]
            selected.append(candidates[0] if len(candidates) == 1 else None)
        kind = row.get("conn_type_id", "").lower()
        if all(selected) and kind in ("covale", "disulf", "metalc", "hydrog"):
            add(
                selected[0],
                selected[1],
                normalized(row.get("pdbx_value_order")) or "unknown",
                "source_" + kind,
                "Explicit retained mmCIF connection: " + str(row.get("id")),
                raw_hash,
            )
        else:
            held.append(
                {"connection": row.get("id"), "reason": "Unresolved typed connection endpoints"}
            )
    water_atoms = sum(a["resname"] in ("HOH", "WAT", "DOD") for a in atoms)
    return bonds, {
        "holds": held,
        "water_atoms_present": water_atoms,
        "hydrogen_atoms_present": sum(a["element"] == "H" for a in atoms),
        "water_scope": (
            f"{water_atoms} water atoms retained; no solvent or missing hydrogen reconstruction"
            if water_atoms
            else "Retained view omits water; no solvent or missing hydrogen reconstruction"
        ),
        "native_context": (
            "First retained model; source chains are not assigned a native partner assembly"
        ),
    }


def build(config: dict) -> dict:
    """Write an immutable optional local atom payload without changing original runtime inputs."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    if config.get("schema") != "marse-atomic-import-config-v1":
        raise ValueError("Unsupported atomic import configuration")
    runtime, output = Path(config["runtime"]), Path(config["output"])
    if output.exists() or output.is_symlink():
        raise ValueError("Atomic output already exists; preserve it and use a successor")
    map_raw = (runtime / "data/map.json").read_bytes()
    meta = json.loads((runtime / "SOURCE_MANIFEST.json").read_bytes())
    ledger_raw = (runtime / "residues.json").read_bytes()
    ledger_meta = json.loads((runtime / "RESIDUE_MANIFEST.json").read_bytes())
    if (
        digest(map_raw) != meta["sha256"]
        or digest(ledger_raw) != ledger_meta["sha256"]
        or ledger_meta["source_snapshot_sha256"] != digest(map_raw)
    ):
        raise ValueError("Current runtime source/residue pins differ")
    data, ledger = json.loads(map_raw), json.loads(ledger_raw)
    fallback_raw = checked(config["template_path"], config["template_sha256"])
    fallback = MMCIF2Dict(io.StringIO(fallback_raw.decode()))
    files, entries, holds = {}, {}, []
    for record in config["models"]:
        ident = record["protein_id"]
        try:
            shape = data["shape_traces"][ident]
            if record["view_sha256"] != shape["source_view_sha256"]:
                raise ValueError("Atomic model does not match the current exact fold")
            view = checked(record["view_path"], record["view_sha256"])
            raw = checked(record["raw_path"], record["raw_sha256"])
            atoms, omitted = parse_pdb(view)
            source_format = record.get("source_format", "mmcif")
            qualification = None
            if source_format == "pdb_prediction":
                if raw != view:
                    raise ValueError(
                        "Retained prediction source must be the exact displayed PDB bytes"
                    )
                protein = next(p for p in data["proteins"] if p["id"] == ident)
                qualification = qualify_pdb_prediction(
                    raw, atoms, shape, ledger["entries"][ident], protein, record["context"]
                )
                doc = {}
            elif source_format == "mmcif":
                doc = MMCIF2Dict(io.StringIO(raw.decode()))
                entry_id = doc.get("_entry.id", [""])[0]
                structure_id = record["context"]["id"]
                predicted_fragment = (
                    record["context"].get("kind") == "prediction"
                    and structure_id.startswith("AF-")
                    and entry_id.startswith(structure_id + "-F")
                    and entry_id[len(structure_id) + 2 :].isdigit()
                    and int(entry_id[len(structure_id) + 2 :]) > 0
                )
                if entry_id != structure_id and not predicted_fragment:
                    raise ValueError("Atomic source entry differs from the retained model")
                qualify_atoms(atoms, doc, shape, ledger["entries"][ident])
            else:
                raise ValueError("Unsupported retained atomic source format")
            bonds, limits = build_topology(atoms, doc, fallback, digest(raw), digest(fallback_raw))
            if source_format == "pdb_prediction":
                limits["source_connections"] = (
                    "PDB prediction has no mmCIF connection categories; "
                    "no disulfides or partners inferred"
                )
                for line in raw.decode().splitlines():
                    if line.startswith(("CONECT", "SSBOND", "LINK  ")):
                        limits["holds"].append(
                            {
                                "reason": (
                                    "PDB connection record not admitted by prediction adapter"
                                ),
                                "connection": line[:6].strip(),
                            }
                        )
            payload = {
                "schema": "marse-atomic-context-v1",
                "protein_id": ident,
                "source_view_sha256": digest(view),
                "raw_source_sha256": digest(raw),
                "source_format": source_format,
                "qualification": qualification,
                "context": record["context"],
                "atoms": atoms,
                "bonds": bonds,
                "selection": {"model": "first", "conformer": "primary/alternate A", **omitted},
                "limitations": limits,
                "scientific_execution": False,
            }
            encoded = (json.dumps(payload, separators=(",", ":"), allow_nan=False) + "\n").encode()
            name = digest(encoded) + ".json"
            files[name] = encoded
            pdb_name = digest(view) + ".pdb"
            files[pdb_name] = view
            entries[ident] = {
                "file": name,
                "sha256": digest(encoded),
                "bytes": len(encoded),
                "pdb_file": pdb_name,
                "source_view_sha256": digest(view),
                "atoms": len(atoms),
                "bonds": len(bonds),
            }
        except (ValueError, KeyError, TypeError) as error:
            holds.append({"protein_id": ident, "reason": str(error)})
    for ident in config.get("requested_proteins", []):
        if ident not in entries and not any(h["protein_id"] == ident for h in holds):
            holds.append(
                {
                    "protein_id": ident,
                    "reason": (
                        "No retained matching full-coordinate structure; "
                        "no atoms or bonds generated"
                    ),
                }
            )
    index = {
        "schema": "marse-atomic-index-v1",
        "source_snapshot_sha256": digest(map_raw),
        "residue_ledger_sha256": digest(ledger_raw),
        "entries": entries,
        "holds": holds,
        "scientific_execution": False,
    }
    index_raw = (json.dumps(index, indent=2) + "\n").encode()
    files["index.json"] = index_raw
    files["IMPORT.json"] = (
        json.dumps(
            {
                "schema": "marse-atomic-import-v1",
                "index_sha256": digest(index_raw),
                "files": {k: digest(v) for k, v in files.items()},
            },
            indent=2,
        )
        + "\n"
    ).encode()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".atomic-import-", dir=output.parent) as temp:
        staging = Path(temp)
        for name, raw in files.items():
            (staging / name).write_bytes(raw)
        staging.rename(output)
    return index


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    result = build(json.loads(args.config.read_bytes()))
    print(
        json.dumps(
            {
                "imported_models": len(result["entries"]),
                "holds": result["holds"],
                "scientific_execution": False,
            }
        )
    )
