"""Import JSON from a retained local MARSE presentation; never execute supplied code."""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
from pathlib import Path, PurePosixPath
from urllib.parse import parse_qs, urlsplit

PAYLOADS = (
    "SOURCE_MANIFEST.json",
    "RESIDUE_MANIFEST.json",
    "residues.json",
    "data/map.json",
    "protein-relation-map.json",
    "RELATION_BUILD.json",
)
MAX_BYTES = 64 * 1024 * 1024


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def checked_source(root: Path, name: str) -> Path:
    """Names in a supplied manifest cannot escape its directory, including via symlinks."""
    if not isinstance(name, str):
        raise ValueError("Invalid retained asset name")
    relative = PurePosixPath(name)
    if relative.is_absolute() or ".." in relative.parts or "\\" in name:
        raise ValueError("Unsafe retained asset name")
    path = (root / name).resolve(strict=True)
    if not path.is_relative_to(root) or not path.is_file() or path.stat().st_size > MAX_BYTES:
        raise ValueError("Retained asset is outside the build or exceeds the size limit")
    return path


def read_build(source: Path) -> dict[str, bytes]:
    root = source.resolve(strict=True)
    if root.name == "studio":
        root = root.parent
    raw = checked_source(root, "BUILD_MANIFEST.json").read_bytes()
    manifest = json.loads(raw)
    files = manifest.get("files", {})
    if (
        manifest.get("schema") != "marse-protein-relations-build-v1"
        or manifest.get("scientific_execution") is not False
        or not isinstance(files, dict)
        or not 1 <= len(files) <= 256
    ):
        raise ValueError("Unsupported retained presentation manifest")
    # Verify all retained files before selecting JSON; supplied scripts are never installed.
    total = 0
    selected = {}
    payload_names = {"studio/" + name: name for name in PAYLOADS}
    for name, expected in files.items():
        contents = checked_source(root, name).read_bytes()
        total += len(contents)
        if len(contents) > MAX_BYTES or total > 512 * 1024 * 1024 or digest(contents) != expected:
            raise ValueError("Retained build checksum or total size differs")
        if name in payload_names:
            selected[payload_names[name]] = contents
    for name in PAYLOADS:
        if name not in selected:
            raise ValueError("Missing retained JSON payload: " + name)
    trace = json.loads(selected["RELATION_BUILD.json"])
    if trace.get("schema") != "marse-molecular-trace-v1":
        raise ValueError("Unsupported retained trace")
    for key in ("build_id", "study_revision_sha256", "packet_sha256"):
        if not manifest.get(key) or trace.get(key) != manifest[key]:
            raise ValueError("Retained manifest/trace identity differs: " + key)
    link = urlsplit(trace["workbench_url"])
    query = parse_qs(link.query)
    if (
        link.scheme != "http"
        or link.hostname not in {"127.0.0.1", "localhost", "::1"}
        or link.username is not None
        or link.password is not None
        or query.get("study") != [manifest["study_id"]]
        or query.get("revision") != [manifest["study_revision_sha256"]]
    ):
        raise ValueError("Unsafe or mismatched exact Workbench return link")
    packet = json.loads(selected["protein-relation-map.json"])
    source_map = json.loads(selected["data/map.json"])
    proteins = {p["id"]: p for p in packet.get("proteins", [])}
    models = {
        key: p["structure_sha256"] for key, p in proteins.items() if p.get("structure_sha256")
    }
    if manifest.get("model_hashes") != models or packet.get("relationships") != source_map.get(
        "relationships"
    ):
        raise ValueError("Retained model or literal-relationship identity differs")
    if source_map.get("subtractive") != packet.get("subtractive_results"):
        raise ValueError("Retained result metadata differs from its packet")
    for native in source_map.get("proteins", []):
        entry = proteins.get(native["id"], {})
        if native.get("subtractive_record") != entry.get("subtractive_record"):
            raise ValueError("Retained result values differ from their packet")
    for key, shape in source_map.get("shape_traces", {}).items():
        if key not in models or shape.get("source_view_sha256") != models[key]:
            raise ValueError("Retained shape/model binding differs")
    source_meta = json.loads(selected["SOURCE_MANIFEST.json"])
    residue_meta = json.loads(selected["RESIDUE_MANIFEST.json"])
    if (
        packet.get("schema") != "protein-relation-map-v1"
        or digest(selected["protein-relation-map.json"]) != manifest["packet_sha256"]
        or digest(selected["data/map.json"]) != source_meta["sha256"]
        or packet["source_snapshot_sha256"] != source_meta["sha256"]
        or digest(selected["residues.json"]) != residue_meta["sha256"]
        or packet["residue_ledger_sha256"] != residue_meta["sha256"]
        or residue_meta["source_snapshot_sha256"] != source_meta["sha256"]
    ):
        raise ValueError("Retained evidence payload bindings differ")
    selected["BUILD_MANIFEST.json"] = raw
    selected["IMPORT.json"] = (
        json.dumps(
            {
                "schema": "marse-local-relations-import-v1",
                "build_id": manifest["build_id"],
                "packet_sha256": manifest["packet_sha256"],
                "retained_manifest_sha256": digest(raw),
                "files": {name: digest(value) for name, value in selected.items()},
                "qualification": (
                    "Retained presentation import; scientific qualification stays upstream."
                ),
            },
            indent=2,
            sort_keys=True,
        ).encode()
        + b"\n"
    )
    return selected


def install(source: Path, destination: Path) -> dict[str, bytes]:
    selected = read_build(source)
    if destination.exists():
        actual = {str(p.relative_to(destination)) for p in destination.rglob("*") if p.is_file()}
        if actual != set(selected) or any(
            not (destination / name).is_file() or (destination / name).read_bytes() != raw
            for name, raw in selected.items()
        ):
            raise ValueError("Existing runtime differs; preserve it and choose a fresh destination")
        return selected
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".relations-import-", dir=destination.parent) as work:
        stage = Path(work) / "payload"
        stage.mkdir()
        for name, raw in selected.items():
            path = stage / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(raw)
        # Rename refuses a nonempty destination created concurrently; no old data is removed.
        stage.rename(destination)
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from", dest="source", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "runtime")
    args = parser.parse_args()
    try:
        selected = install(args.source, args.output)
    except (ValueError, KeyError, OSError, json.JSONDecodeError) as exc:
        raise SystemExit("Local presentation import refused: " + str(exc)) from exc
    print(f"Checked {len(selected)} local JSON payloads; no supplied code or network access.")


if __name__ == "__main__":
    main()
