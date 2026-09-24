"""Tests for reproducible ecosystem state-space ensembles."""

import json

import pytest

from marse.ensemble import RunCatalog, ScenarioBatch, run_batch


def base() -> dict:
    return {
        "experiment_id": "ensemble",
        "width": 4,
        "height": 4,
        "cell_size_um": 10,
        "duration_h": 0.1,
        "timestep_h": 0.01,
        "seed": 1,
        "nutrients": [{"name": "food", "initial": 1, "diffusivity": 0}],
        "species": [
            {
                "name": "one",
                "initial_biomass": 0.1,
                "maximum_growth_per_h": 0.2,
                "half_saturation": [0.2],
                "yield_per_nutrient": [1],
            }
        ],
    }


def test_grid_is_content_addressed_and_deduplicated():
    batch = ScenarioBatch.grid(base(), {"nutrients.0.initial": [1, 2, 1]})
    assert len(batch.scenarios) == 2
    assert batch.scenarios[0].scenario_id.startswith("SCN-")


def test_batch_execution_returns_compact_signatures_and_is_reproducible(tmp_path):
    batch = ScenarioBatch.grid(base(), {"seed": [1, 2]})
    first = run_batch(batch)
    second = run_batch(batch, workers=2)
    assert first.records == second.records
    assert all(record.status == "completed" for record in first.records)
    assert "biomass:one" in first.records[0].signature
    path = RunCatalog(first.records).write_json(tmp_path / "catalog.json")
    assert json.loads(path.read_text())["runs"][0]["scenario_id"] == batch.scenarios[0].scenario_id


def test_invalid_worker_count_and_failed_scenario_are_explicit():
    with pytest.raises(ValueError, match="workers"):
        run_batch(ScenarioBatch.grid(base(), {}), workers=0)
    invalid = ScenarioBatch.grid(base(), {"width": [1]})
    record = run_batch(invalid).records[0]
    assert record.status == "failed"
    assert record.error is not None
