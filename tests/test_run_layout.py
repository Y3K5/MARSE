"""Every engine leaves the same directory records, and old manifests still replay.

``summary.json`` and ``report.html`` sit beside the scientific files. They are
not part of the manifest, so the configuration checksum is the one the schema
has always produced.
"""

from __future__ import annotations

import json
from pathlib import Path

from marse.cli import main
from marse.core.provenance import Manifest


def _batch(path: Path) -> None:
    path.write_text(
        json.dumps(
            {
                "experiment_id": "layout-batch",
                "seed": 1,
                "duration_h": 0.01,
                "timestep_h": 0.01,
                "checkpoint_interval_h": 0.01,
                "environment": {"temperature_c": 37.0, "ph": 7.0},
                "substrate": {"name": "glucose", "initial_mm": 10.0},
                "organisms": [
                    {
                        "name": "solo",
                        "initial_biomass_g_per_l": 0.01,
                        "mu_opt_per_h": 0.8,
                        "k_s_mm": 0.4,
                        "yield_g_per_mmol": 0.08,
                    }
                ],
            }
        ),
        encoding="utf-8",
    )


def test_a_run_directory_has_the_shared_records_and_an_unchanged_manifest(tmp_path: Path):
    experiment = tmp_path / "batch.json"
    _batch(experiment)
    output = tmp_path / "store" / "batch"
    assert main(["run", str(experiment), "-o", str(output)]) == 0
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert "layout_version" not in manifest
    assert "summary.json" not in json.dumps(manifest)
    assert summary["layout_version"] == 1
    assert summary["engine"] == "legacy_batch"
    assert summary["kind"] == "batch"
    assert summary["config_sha256"] == manifest["config_sha256"]
    assert summary["unverified"] is False
    assert (output / "trajectory.csv").is_file()
    assert (output / "report.html").is_file()
    assert "layout-batch" in (output / "report.html").read_text(encoding="utf-8")
    index = json.loads((output.parent / "index.json").read_text(encoding="utf-8"))
    assert index["runs"][0]["directory"] == "batch"
    assert index["runs"][0]["config_sha256"] == manifest["config_sha256"]
    assert str(output) not in (output.parent / "index.json").read_text(encoding="utf-8")
    assert main(["replay", str(output / "manifest.json")]) == 0
    assert main(["ls", str(output.parent)]) == 0
    other = tmp_path / "store" / "again"
    assert main(["run", str(experiment), "-o", str(other)]) == 0
    assert main(["compare", str(output), str(other)]) == 0


def test_compare_reports_a_different_configuration(tmp_path: Path):
    first = tmp_path / "a.json"
    second = tmp_path / "b.json"
    _batch(first)
    _batch(second)
    document = json.loads(second.read_text(encoding="utf-8"))
    document["seed"] = 2
    second.write_text(json.dumps(document), encoding="utf-8")
    assert main(["run", str(first), "-o", str(tmp_path / "one")]) == 0
    assert main(["run", str(second), "-o", str(tmp_path / "two")]) == 0
    assert main(["compare", str(tmp_path / "one" / "manifest.json"), str(tmp_path / "two")]) == 1
    # Reading the manifest the way a previous release wrote it still works.
    Manifest.read(tmp_path / "one" / "manifest.json")
