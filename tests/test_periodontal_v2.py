"""The periodontal study on the version-2 engine.

The file is a column of three obligate anaerobes sharing alanine. It does not
reproduce the version-1 cell map, and it does not claim a calibrated pocket.
The version-1 experiment is still there, and still unverified.
"""

from __future__ import annotations

import json
from pathlib import Path

from marse.cli import main
from marse.engines import registry

EXAMPLE = Path("examples/networks/periodontal_pathogen_v2.json")
LEGACY = Path("examples/experiments/periodontal_pathogen_biofilm.json")


def test_the_v2_study_is_three_anaerobes_on_the_spatial_engine():
    engine, plan = registry.load(EXAMPLE)
    assert engine.name == "reactive_transport"
    roles = {
        component.name: component.oxygen_role
        for component in plan.network.components
        if component.oxygen_role is not None
    }
    assert roles == {
        "t_denticola": "obligate_anaerobe",
        "t_forsythia": "obligate_anaerobe",
        "p_gingivalis": "obligate_anaerobe",
    }
    assert "oxygen" not in {component.name for component in plan.network.components}
    for process in plan.network.processes:
        names = set(process.given) | set(process.balanced_by)
        if process.rate is not None:
            names.update(factor.component for factor in process.rate.factors)
        assert "oxygen" not in names


def test_the_v1_study_is_still_the_unverified_ecosystem():
    """``marse ecosystem`` still loads the version-1 file. ``marse run`` does not infer it."""
    from marse.engines.legacy_ecosystem import legacy_ecosystem_v1

    raw = json.loads(LEGACY.read_text(encoding="utf-8"))
    assert "schema_version" not in raw
    assert registry.infer(raw).name == "legacy_batch"
    plan = legacy_ecosystem_v1.load(LEGACY)
    assert plan.experiment_id == "periodontal-pathogen-biofilm"
    assert legacy_ecosystem_v1.unverified


def test_the_v2_column_runs_and_replays(tmp_path: Path):
    output = tmp_path / "v2"
    assert main(["run", str(EXAMPLE), "-o", str(output)]) == 0
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["kind"] == "reactive_transport"
    balance = manifest["outputs"]["balance"]
    for quantity in balance.values():
        assert quantity["largest_relative_residual"] <= 1e-12
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert summary["engine"] == "reactive_transport"
    assert summary["experiment_id"] == "periodontal-pathogen-v2"
    assert main(["replay", str(output / "manifest.json"), "-o", str(tmp_path / "replay")]) == 0
