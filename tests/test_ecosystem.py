"""Tests for the deterministic 2D multi-species ecosystem foundation."""

from pathlib import Path

import numpy as np
import pytest

from marse.ecosystem import (
    EcosystemConfig,
    EcosystemError,
    NutrientConfig,
    SpeciesConfig,
    load_experiment,
    run,
    write_viewer,
)


def config(**overrides) -> EcosystemConfig:
    values = {
        "experiment_id": "ecosystem-test",
        "width": 8,
        "height": 6,
        "cell_size_um": 10.0,
        "duration_h": 1.0,
        "timestep_h": 0.01,
        "seed": 42,
        "nutrients": (
            NutrientConfig("oxygen", 1.0, 10.0),
            NutrientConfig("carbon", 2.0, 5.0),
        ),
        "species": (
            SpeciesConfig("aerobe", 0.1, 0.8, (0.2, 0.4), (1.0, 1.0)),
            SpeciesConfig("competitor", 0.1, 0.4, (0.3, 0.2), (1.0, 1.0)),
        ),
    }
    return EcosystemConfig(**(values | overrides))


def test_species_grow_and_nutrients_are_consumed():
    result = run(config())
    assert result.final_state.time_h == pytest.approx(1.0)
    assert result.final_state.biomass.shape == (2, 6, 8)
    assert np.all(result.final_state.biomass > 0.1)
    assert np.all(result.final_state.nutrients < np.array([1.0, 2.0])[:, None, None])


def test_same_seed_reproduces_all_frames():
    first, second = run(config(seed=9)), run(config(seed=9))
    for a, b in zip(first.frames, second.frames, strict=True):
        np.testing.assert_array_equal(a.biomass, b.biomass)
        np.testing.assert_array_equal(a.nutrients, b.nutrients)
        np.testing.assert_array_equal(a.mutations, b.mutations)


def test_mutations_are_seeded_and_recorded():
    result = run(
        config(
            duration_h=0.2,
            timestep_h=0.01,
            mutation_interval_h=0.1,
            species=(SpeciesConfig("mutant", 0.1, 0.2, (0.2,), (1.0,), 1.0, 2.0),),
            nutrients=(NutrientConfig("food", 1.0, 0.0),),
        )
    )
    assert result.final_state.mutations.sum() > 0
    assert result.frames[-1].mutations.shape == (1, 6, 8)


def test_diffusion_cfl_is_checked():
    with pytest.raises(EcosystemError, match="unstable"):
        run(config(timestep_h=3.0, duration_h=3.0))


def test_frames_can_be_exported_for_a_visualizer(tmp_path: Path):
    path = run(config()).write_frames(tmp_path / "frames.json")
    text = path.read_text()
    assert '"species"' in text
    assert '"oxygen"' in text
    assert '"time_h":0.0' in text


def test_json_experiment_loader_and_viewer(tmp_path: Path):
    source = tmp_path / "experiment.json"
    source.write_text(
        '{"experiment_id":"json","width":4,"height":4,"cell_size_um":10,'
        '"duration_h":0.1,"timestep_h":0.01,"seed":1,'
        '"nutrients":[{"name":"food","initial":1,"diffusivity":0}],'
        '"species":[{"name":"one","initial_biomass":0.1,"maximum_growth_per_h":0.2,'
        '"half_saturation":[0.2],"yield_per_nutrient":[1]}]}'
    )
    result = run(load_experiment(source))
    viewer = write_viewer(result, tmp_path / "viewer.html")
    assert viewer.is_file()
    assert (tmp_path / "frames.json").is_file()
    assert "canvas" in viewer.read_text()
