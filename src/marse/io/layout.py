"""One run directory for every engine, plus an index of the directories beside it.

Every engine still writes the scientific files it always has (a trajectory, a
profile, totals, frames, a viewer). On top of those, each run directory gains
the same three records:

- ``manifest.json``, unchanged, so an old manifest still replays;
- ``summary.json``, which names the engine and lists the files;
- ``report.html``, a short page that states what the run is and is not.

The manifest is the checksummed record. ``summary.json`` and ``report.html``
are not part of it, so adding them does not change ``config_sha256`` or the
recorded outputs. The parent directory's ``index.json`` lists those runs by
folder name. It stores no absolute path.
"""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

from marse.core.provenance import Manifest

LAYOUT_VERSION = 1
INDEX_VERSION = 1
_SUMMARY = "summary.json"
_REPORT = "report.html"
_INDEX = "index.json"


def publish_run(output_dir: Path, manifest: Manifest) -> None:
    """Write ``summary.json``, ``report.html``, and update the parent index."""
    output_dir.mkdir(parents=True, exist_ok=True)
    engine = _engine_name(manifest.kind)
    unverified = _engine_unverified(manifest.kind)
    report = output_dir / _REPORT
    report.write_text(_report_html(manifest, engine, unverified), encoding="utf-8")
    files = sorted(
        path.relative_to(output_dir).as_posix()
        for path in output_dir.rglob("*")
        if path.is_file() and path.name != _SUMMARY
    )
    files.append(_SUMMARY)
    summary = {
        "layout_version": LAYOUT_VERSION,
        "run_id": manifest.run_id,
        "experiment_id": manifest.experiment_id,
        "kind": manifest.kind,
        "engine": engine,
        "unverified": unverified,
        "config_sha256": manifest.config_sha256,
        "started_utc": manifest.started_utc,
        "finished_utc": manifest.finished_utc,
        "files": files,
    }
    (output_dir / _SUMMARY).write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _update_index(output_dir, summary)


def list_runs(directory: Path) -> list[dict[str, Any]]:
    """Runs under ``directory``, from manifests, so a directory without an index still lists."""
    directory = directory.resolve()
    found: list[dict[str, Any]] = []
    if (directory / "manifest.json").is_file():
        manifests = [directory / "manifest.json"]
    else:
        manifests = []
        for path in _walk_manifests(directory, depth=3):
            manifests.append(path)
    for path in manifests:
        manifest = Manifest.read(path)
        run_dir = path.parent
        try:
            folder = run_dir.relative_to(directory).as_posix()
        except ValueError:
            folder = run_dir.name
        if folder == ".":
            folder = run_dir.name
        found.append(
            {
                "directory": folder,
                "run_id": manifest.run_id,
                "experiment_id": manifest.experiment_id,
                "kind": manifest.kind,
                "engine": _engine_name(manifest.kind),
                "config_sha256": manifest.config_sha256,
                "finished_utc": manifest.finished_utc,
            }
        )
    return sorted(found, key=lambda row: row["directory"])


def read_manifest(path: Path) -> Manifest:
    """A manifest file, or ``manifest.json`` inside a run directory."""
    if path.is_dir():
        path = path / "manifest.json"
    return Manifest.read(path)


def _engine_name(kind: str) -> str:
    from marse.engines.registry import registry

    try:
        return registry.for_kind(kind).name
    except ValueError:
        return ""


def _engine_unverified(kind: str) -> bool:
    from marse.engines.registry import registry

    try:
        return bool(registry.for_kind(kind).unverified)
    except ValueError:
        return False


def _update_index(output_dir: Path, summary: dict[str, Any]) -> None:
    parent = output_dir.resolve().parent
    index_path = parent / _INDEX
    if index_path.is_file():
        try:
            raw = json.loads(index_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            raw = {}
    else:
        raw = {}
    runs = [row for row in raw.get("runs", []) if row.get("directory") != output_dir.name]
    runs.append(
        {
            "directory": output_dir.name,
            "run_id": summary["run_id"],
            "experiment_id": summary["experiment_id"],
            "kind": summary["kind"],
            "engine": summary["engine"],
            "unverified": summary["unverified"],
            "config_sha256": summary["config_sha256"],
            "finished_utc": summary["finished_utc"],
        }
    )
    document = {
        "index_version": INDEX_VERSION,
        "runs": sorted(runs, key=lambda row: row["directory"]),
    }
    index_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _walk_manifests(directory: Path, depth: int) -> list[Path]:
    found: list[Path] = []

    def walk(current: Path, remaining: int) -> None:
        if remaining < 0 or not current.is_dir():
            return
        manifest = current / "manifest.json"
        if manifest.is_file() and current != directory:
            found.append(manifest)
            return
        for child in sorted(current.iterdir()):
            if child.is_dir() and child.name not in {"frames", "vtk", "__pycache__"}:
                walk(child, remaining - 1)

    walk(directory, depth)
    return found


def _report_html(manifest: Manifest, engine: str, unverified: bool) -> str:
    title = html.escape(manifest.experiment_id)
    kind = html.escape(manifest.kind)
    engine_name = html.escape(engine or "unknown")
    checksum = html.escape(manifest.config_sha256)
    note = (
        "This run used the unverified 2-D ecosystem engine. "
        "docs/validation.md records known defects. The numbers support no quantitative conclusion."
        if unverified
        else "Replay proves this run can be repeated. It does not prove the model is right."
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
body {{ font-family: sans-serif; max-width: 40rem; margin: 2rem auto; line-height: 1.4; }}
dt {{ font-weight: bold; }}
</style>
</head>
<body>
<h1>{title}</h1>
<dl>
<dt>Run</dt><dd>{html.escape(manifest.run_id)}</dd>
<dt>Kind</dt><dd>{kind}</dd>
<dt>Engine</dt><dd>{engine_name}</dd>
<dt>Config checksum</dt><dd><code>{checksum}</code></dd>
</dl>
<p>{note}</p>
</body>
</html>
"""
