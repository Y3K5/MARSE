"""Create an explicitly invented, local software fixture, never a biological example."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def raw(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode() + b"\n"


def sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def make_fixture(root: Path) -> dict:
    if root.exists():
        raise ValueError("Fixture output already exists")
    species = [
        {"id": name, "name": "Synthetic " + name, "short": name.title(), "form": "test marker"}
        for name in ("alpha", "beta", "gamma")
    ]
    proteins, shapes, entries, packet_proteins, placements = [], {}, {}, [], []
    rows = [
        ("review-a", "Review A", "alpha", "A", True, True),
        ("review-hold", "Probe review", "alpha", "G", True, True),
        ("review-gap", "Coordinate gap", "alpha", "V", False, True),
        ("source-a", "Probe", "alpha", "A", True, False),
        ("source-b", "Source B", "beta", "L", True, False),
        ("source-c", "Source C", "gamma", "S", True, False),
    ]
    residue_names = {"A": "ALA", "G": "GLY", "V": "VAL", "L": "LEU", "S": "SER"}
    for ident, name, taxon, letter, folded, result in rows:
        sequence = letter * 8
        sequence_hash = sha(sequence.encode())
        model_hash = sha(("invented-model-" + ident).encode()) if folded else None
        protein = {
            "id": ident,
            "name": name,
            "accession": "SYN-" + ident,
            "species": taxon,
            "kind": "protein",
            "category": "Synthetic fixture",
            "role": "Software test only",
            "localization": "Unresolved",
            "caveat": "Invented, not biological evidence",
            "architecture": [],
            "source_ids": ["synthetic"],
            "reference": {
                "sequence": sequence,
                "length": 8,
                "organism": {"scientificName": "Invented " + taxon},
            },
        }
        if result:
            record = {
                "key": "alpha-ref|" + ident,
                "id": "SYN-" + ident,
                "name": name,
                "alias": name,
                "species": "alpha-ref",
                "strain": "Invented reference",
                "length": 8,
                "sequence_sha256": sequence_hash,
                "primary_group": "unresolved",
                "source_rank": 1,
                "score": None,
                "review_label": "Synthetic review only",
                "holds": ["Software fixture; never adopt as a target"],
            }
            protein.update(
                subtractive_record=record,
                catalogue_key=record["key"],
                retained_evidence={
                    "structure": {
                        "status": "exact_prediction" if folded else "missing",
                        "label": "Synthetic coordinates" if folded else "Gap",
                    },
                    "location": {"label": "Unresolved"},
                },
            )
        proteins.append(protein)
        placements.append({"protein_id": ident, "zone": "unresolved", "note": "Display only"})
        packet_proteins.append(
            {
                "id": ident,
                "accession": protein["accession"],
                "sequence_sha256": sequence_hash,
                "structure_sha256": model_hash,
                "strain_reference": protein["reference"]["organism"]["scientificName"],
            }
        )
        if result:
            packet_proteins[-1]["subtractive_record"] = protein["subtractive_record"]
        if folded:
            coordinates = [[i + 1, i * 2.0, math.sin(i) * 2, math.cos(i) * 2] for i in range(8)]
            shapes[ident] = {
                "id": "SYNTHETIC-" + ident,
                "kind": "prediction",
                "state": "Invented test input",
                "scope": "Synthetic coordinates, not a predicted or measured biological fold",
                "source_view_sha256": model_hash,
                "chains": [{"chain": "A", "ca": coordinates}],
            }
            entries[ident] = {
                "source_view_sha256": model_hash,
                "chains": [
                    {
                        "chain": "A",
                        "residues": [[i + 1, residue_names[letter], ""] for i in range(8)],
                    }
                ],
            }
    relationships = [
        {
            "id": "pair-1",
            "a": "source-a",
            "b": "source-b",
            "relation": "synthetic pair",
            "level": "direct assay",
            "source_ids": ["synthetic"],
            "direction": "undirected",
            "claim": "Invented pair for software testing; not measured evidence",
            "limit": "No biological inference",
        },
        {
            "id": "candidate-1",
            "a": "source-a",
            "b": "gamma-cell",
            "relation": "synthetic candidate",
            "level": "candidate / unresolved identity",
            "source_ids": ["synthetic"],
            "direction": "directed",
            "claim": "Invented unresolved partner",
            "limit": "No partner assigned",
        },
    ]
    data = {
        "schema": "grouped-organism-workspace-v1",
        "native_placements_admitted": 0,
        "candidate_interaction_transfer": False,
        "presentation": {"testing_only": True},
        "species": species,
        "proteins": proteins,
        "shape_traces": shapes,
        "placements": placements,
        "contexts": [{"id": "gamma-cell", "name": "Unresolved synthetic cell", "species": "gamma"}],
        "systems": [],
        "relationships": relationships,
        "pairs": {},
        "workspace_bridge": {},
        "study_contexts": {},
        "sources": [
            {
                "id": "synthetic",
                "title": "Invented fixture, no paper",
                "url": "https://example.invalid/software-fixture",
                "year": "Not applicable",
                "evidence_note": "Software test only",
            }
        ],
        "subtractive": {
            "schema": "retained-subtractive-marse-v1",
            "record_count": 3,
            "saved_form_count": 2,
            "exact_fold_count": 2,
            "source_context_ids": ["source-a", "source-b", "source-c"],
            "shortlist_keys": [
                p["subtractive_record"]["key"] for p in proteins if "subtractive_record" in p
            ],
            "groups": [{"id": "unresolved", "label": "Unresolved"}],
            "species": [
                {
                    "id": "alpha-ref",
                    "label": "Synthetic reference",
                    "rank_origin": "Invented",
                    "rank_identity_limit": "Software fixture, no scientific execution",
                }
            ],
            "reference_links": [
                {
                    "id": "reference-1",
                    "a": "review-a",
                    "b": "source-a",
                    "kind": "exact_reference_identity",
                    "sequence_sha256": sha(b"AAAAAAAA"),
                }
            ],
        },
    }
    source_bytes = raw(data)
    ledger = {
        "schema": "marse-source-residues-v1",
        "source_snapshot_sha256": sha(source_bytes),
        "entries": entries,
    }
    ledger_bytes = raw(ledger)
    annotation = {
        "id": "synthetic-construct",
        "finding_id": "synthetic-finding",
        "relationship_id": "pair-1",
        "endpoint_id": "source-a",
        "model_protein_id": "source-a",
        "category": "assayed_construct",
        "status": "supported_location",
        "default_visible": True,
        "exact_contact": False,
        "counterpart_assignment": None,
        "label": "Synthetic test span, not a biological construct",
        "limitations": "Invented fixture only; no contact site or partner coordinates",
        "identity": {
            "accession": "SYN-source-a",
            "structure_sha256": shapes["source-a"]["source_view_sha256"],
            "sequence_sha256": sha(b"AAAAAAAA"),
        },
        "numbering": {"paper": "Synthetic numbering"},
        "correspondence": [
            {"paper": i, "canonical": i, "author": i, "chain": "A", "insertion": "", "name": "ALA"}
            for i in (2, 3)
        ],
    }
    packet = {
        "schema": "protein-relation-map-v1",
        "subtractive_results": data["subtractive"],
        "source_snapshot_sha256": sha(source_bytes),
        "residue_ledger_sha256": sha(ledger_bytes),
        "relationships": relationships,
        "proteins": packet_proteins,
        "annotations": [annotation],
        "sources": data["sources"],
        "findings": [
            {
                "id": "synthetic-finding",
                "material": "Invented software input; no biological assay",
                "locator": "Synthetic fixture, not a publication",
                "method": {"name": "Software fixture", "version": "fixture-v1"},
                "endpoint_id": "source-a",
                "category": "assayed_construct",
                "source_ids": ["synthetic"],
            }
        ],
        "cross_microbe_audit": [],
    }
    packet_bytes = raw(packet)
    core = {
        "schema": "marse-protein-relations-build-v1",
        "study_id": "synthetic-study",
        "study_revision_sha256": "b" * 64,
        "packet_sha256": sha(packet_bytes),
        "renderer_version": "synthetic-fixture-v1",
        "renderer_sha256": "c" * 64,
        "layout": "cluster",
        "scientific_execution": False,
        "model_hashes": {
            p["id"]: p["structure_sha256"] for p in packet_proteins if p["structure_sha256"]
        },
    }
    identity = sha(raw(core))
    trace = {
        "schema": "marse-molecular-trace-v1",
        "build_id": identity,
        "study_revision_sha256": core["study_revision_sha256"],
        "packet_sha256": core["packet_sha256"],
        "workbench_url": "http://127.0.0.1:8990/molecular-relations/?study=synthetic-study&revision="
        + "b" * 64,
    }
    payloads = {
        "SOURCE_MANIFEST.json": raw({"sha256": sha(source_bytes)}),
        "RESIDUE_MANIFEST.json": raw(
            {
                "sha256": sha(ledger_bytes),
                "source_snapshot_sha256": sha(source_bytes),
                "residue_count": 40,
            }
        ),
        "data/map.json": source_bytes,
        "residues.json": ledger_bytes,
        "protein-relation-map.json": packet_bytes,
        "RELATION_BUILD.json": raw(trace),
    }
    root.mkdir(parents=True)
    for name, value in payloads.items():
        path = root / "studio" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(value)
    manifest = {
        **core,
        "build_id": identity,
        "files": {"studio/" + name: sha(value) for name, value in payloads.items()},
    }
    (root / "BUILD_MANIFEST.json").write_bytes(raw(manifest))
    return manifest


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    make_fixture(parser.parse_args().output)
    print("Created invented software fixture; not biological evidence.")
