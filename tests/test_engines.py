"""Engines are chosen from the experiment file, not from branches in the command.

An explicit ``engine`` field selects one and is removed before the scientific
schema sees the document, so it does not change the manifest. Without it,
``marse run`` recognises a file the way it always has.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from marse.cli import main
from marse.core.config import ConfigError
from marse.engines import registry
from marse.engines.legacy_ecosystem import UNVERIFIED_WARNING

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def test_a_version_2_network_without_a_domain_is_well_mixed():
    engine, plan = registry.load(EXAMPLES / "networks" / "glucose_cross_feeding.json")
    assert engine.name == "well_mixed"
    assert plan.kind == "well_mixed"


def test_a_version_2_network_with_a_domain_is_reactive_transport():
    engine, plan = registry.load(EXAMPLES / "networks" / "surface_biofilm_1d.json")
    assert engine.name == "reactive_transport"
    assert plan.kind == "reactive_transport"


def test_a_batch_experiment_and_a_biofilm_profile_select_their_own_engines():
    batch, _ = registry.load(EXAMPLES / "experiments" / "two_species_batch.json")
    biofilm, _ = registry.load(EXAMPLES / "experiments" / "biofilm_oxygen_profile.json")
    assert batch.name == "legacy_batch"
    assert biofilm.name == "steady_biofilm"


def test_an_explicit_engine_is_stripped_before_the_schema_and_replay_still_matches(tmp_path: Path):
    document = {
        "engine": "legacy_batch",
        "experiment_id": "engine-field",
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
    named = tmp_path / "named.json"
    plain = tmp_path / "plain.json"
    named.write_text(json.dumps(document), encoding="utf-8")
    body = {key: value for key, value in document.items() if key != "engine"}
    plain.write_text(json.dumps(body), encoding="utf-8")
    named_out, plain_out = tmp_path / "named", tmp_path / "plain"
    assert main(["run", str(named), "-o", str(named_out)]) == 0
    assert main(["run", str(plain), "-o", str(plain_out)]) == 0
    named_manifest = json.loads((named_out / "manifest.json").read_text(encoding="utf-8"))
    plain_manifest = json.loads((plain_out / "manifest.json").read_text(encoding="utf-8"))
    assert "engine" not in named_manifest["config"]
    assert named_manifest["config_sha256"] == plain_manifest["config_sha256"]
    assert named_manifest["outputs"] == plain_manifest["outputs"]


def test_the_named_engine_must_agree_with_the_document():
    path = EXAMPLES / "networks" / "glucose_cross_feeding.json"
    raw = json.loads(path.read_text(encoding="utf-8"))
    raw["engine"] = "reactive_transport"
    with pytest.raises(ConfigError, match="needs a domain"):
        registry.get("reactive_transport").validate(
            {key: value for key, value in raw.items() if key != "engine"}
        )


def test_an_unknown_engine_is_refused(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    path = tmp_path / "nope.json"
    path.write_text(json.dumps({"engine": "not-an-engine"}), encoding="utf-8")
    assert main(["run", str(path)]) == 2
    assert "not registered" in capsys.readouterr().out


def test_the_ecosystem_command_warns_that_the_engine_is_unverified(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
):
    document = {
        "experiment_id": "tiny-ecosystem",
        "width": 4,
        "height": 4,
        "cell_size_um": 10.0,
        "duration_h": 0.01,
        "timestep_h": 0.01,
        "seed": 1,
        "nutrients": [{"name": "carbon", "initial": 1.0, "diffusivity": 1.0}],
        "species": [
            {
                "name": "grower",
                "initial_biomass": 0.1,
                "maximum_growth_per_h": 0.5,
                "half_saturation": [0.2],
                "yield_per_nutrient": [1.0],
            }
        ],
    }
    path = tmp_path / "ecosystem.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    assert main(["ecosystem", str(path), "-o", str(tmp_path / "out")]) == 0
    captured = capsys.readouterr()
    assert UNVERIFIED_WARNING in captured.err
    assert "run_id" in captured.out
    assert (tmp_path / "out" / "manifest.json").is_file()
