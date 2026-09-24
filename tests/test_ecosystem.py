"""Tests for the deterministic 2D multi-species ecosystem foundation."""

from pathlib import Path

import numpy as np
import pytest

from marse.additives import AdditiveEffect, apply_effect, hill_response
from marse.ecosystem import (
    AdditiveConfig,
    ConditionConfig,
    EcosystemConfig,
    EcosystemError,
    EcosystemProviders,
    ExplicitTransportProvider,
    NutrientConfig,
    PhenotypeConfig,
    SeedRegion,
    SpeciesConfig,
    load_experiment,
    run,
    write_viewer,
)
from marse.niche import Capability


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
    stats = result.frames[-1].to_dict(("aerobe", "competitor"), ("oxygen", "carbon"), 1.0)
    assert stats["statistics"]["species_total_biomass"]["aerobe"] > 0


def test_same_seed_reproduces_all_frames():
    first, second = run(config(seed=9)), run(config(seed=9))
    for a, b in zip(first.frames, second.frames, strict=True):
        np.testing.assert_array_equal(a.biomass, b.biomass)
        np.testing.assert_array_equal(a.nutrients, b.nutrients)
        np.testing.assert_array_equal(a.mutations, b.mutations)


def test_default_provider_pipeline_is_declared_and_exported():
    result = run(config())
    assert result.provider_versions == EcosystemProviders().versions
    assert result.provider_versions["transport"] == "explicit_transport_2d_v1"
    assert result.provider_versions["reactions"] == "local_reactions_2d_v1"


def test_custom_transport_provider_preserves_provider_metadata():
    result = run(
        config(),
        EcosystemProviders(transport=ExplicitTransportProvider(version="test_transport_v1")),
    )
    assert result.provider_versions["transport"] == "test_transport_v1"
    assert result.final_state.time_h == pytest.approx(1.0)


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


def test_seed_regions_create_localized_initial_colonies():
    result = run(
        config(
            duration_h=0.1,
            species=(
                SpeciesConfig(
                    "localized",
                    0.0,
                    0.1,
                    (0.2, 0.2),
                    (1.0, 1.0),
                    seed_regions=(SeedRegion(0.25, 0.5, 0.2, 0.8),),
                ),
                config().species[1],
            ),
        )
    )
    field = result.frames[0].biomass[0]
    assert field.max() == pytest.approx(0.8)
    assert np.count_nonzero(field) < field.size


def test_fixed_nutrient_boundary_creates_a_gradient():
    result = run(
        config(
            duration_h=0.5,
            timestep_h=0.01,
            nutrients=(
                NutrientConfig(
                    "oxygen",
                    0.0,
                    10.0,
                    boundary_value=1.0,
                    boundary_edges=("top",),
                ),
                NutrientConfig("carbon", 2.0, 0.0),
            ),
            species=(SpeciesConfig("consumer", 0.0, 0.2, (0.2, 0.2), (1.0, 1.0)),),
        )
    )
    oxygen = result.final_state.nutrients[0]
    assert np.all(oxygen[0] == pytest.approx(1.0))
    assert oxygen[-1].mean() < oxygen[0].mean()


def test_json_spatial_fields_are_loaded(tmp_path: Path):
    source = tmp_path / "spatial.json"
    source.write_text(
        '{"experiment_id":"spatial","width":4,"height":4,"cell_size_um":10,'
        '"duration_h":0.1,"timestep_h":0.01,"seed":1,'
        '"nutrients":[{"name":"food","initial":0,"diffusivity":1,'
        '"boundary_value":1,"boundary_edges":["top"]}],'
        '"species":[{"name":"one","initial_biomass":0,"maximum_growth_per_h":0.2,'
        '"half_saturation":[0.2],"yield_per_nutrient":[1],'
        '"seed_regions":[{"x":0.5,"y":0.5,"radius":0.2,"biomass":1}]}]}'
    )
    loaded = load_experiment(source)
    assert loaded.nutrients[0].boundary_edges == ("top",)
    assert loaded.species[0].seed_regions[0].biomass == 1


def test_competition_limits_growth_and_spreading_reaches_neighbors():
    base_species = (
        SpeciesConfig(
            "one",
            0.0,
            1.0,
            (0.2,),
            (1.0,),
            spreading_per_h=5.0,
            seed_regions=(SeedRegion(0.5, 0.5, 0.2, 0.4),),
        ),
        SpeciesConfig(
            "two",
            0.0,
            1.0,
            (0.2,),
            (1.0,),
            seed_regions=(SeedRegion(0.5, 0.5, 0.2, 0.4),),
        ),
    )
    common = dict(
        duration_h=0.2,
        timestep_h=0.001,
        cell_size_um=10.0,
        nutrients=(NutrientConfig("food", 5.0, 0.0),),
        species=base_species,
        carrying_capacity=0.5,
    )
    without = run(config(**common, competition_coefficients=((0.0, 0.0), (0.0, 0.0))))
    with_competition = run(config(**common, competition_coefficients=((0.0, 4.0), (4.0, 0.0))))
    assert with_competition.final_state.biomass.sum() < without.final_state.biomass.sum()
    assert np.count_nonzero(with_competition.final_state.biomass[0] > 0) > np.count_nonzero(
        with_competition.frames[0].biomass[0] > 0
    )


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


def test_species_can_produce_a_metabolite_for_another_species():
    producer = SpeciesConfig(
        "producer",
        0.2,
        0.8,
        (0.2, 0.01),
        (1.0, 1.0),
        production_per_nutrient=(0.0, 2.0),
    )
    consumer = SpeciesConfig(
        "consumer",
        0.1,
        0.5,
        (10.0, 0.02),
        (1.0, 1.0),
    )
    result = run(
        config(
            duration_h=0.5,
            timestep_h=0.001,
            nutrients=(
                NutrientConfig("carbon", 1.0, 0.0),
                NutrientConfig("acetate", 0.1, 1.0),
            ),
            species=(producer, consumer),
        )
    )
    initial_acetate = result.frames[0].nutrients[1].sum()
    final_acetate = result.final_state.nutrients[1].sum()
    assert final_acetate > initial_acetate
    assert result.final_state.biomass[1].sum() > result.frames[0].biomass[1].sum()


def test_negative_production_is_rejected():
    with pytest.raises(EcosystemError, match="production"):
        SpeciesConfig("invalid", 0.1, 0.2, (0.2,), (1.0,), production_per_nutrient=(-1.0,))


def test_chemotaxis_moves_biomass_up_a_declared_signal_gradient():
    result = run(
        config(
            width=6,
            height=6,
            duration_h=0.2,
            timestep_h=0.001,
            nutrients=(NutrientConfig("food", 1.0, 0.0),),
            conditions=(
                ConditionConfig(
                    "signal",
                    0.0,
                    0.0,
                    boundary_value=1.0,
                    boundary_edges=("top",),
                ),
            ),
            species=(
                SpeciesConfig(
                    "motile",
                    0.0,
                    0.0,
                    (0.2,),
                    (1.0,),
                    chemotaxis_field="signal",
                    chemotaxis_sensitivity=100.0,
                    seed_regions=(SeedRegion(0.5, 0.2, 0.12, 1.0),),
                ),
            ),
        )
    )
    initial_center = np.average(
        np.indices(result.frames[0].biomass[0].shape)[0],
        weights=result.frames[0].biomass[0],
    )
    final_center = np.average(
        np.indices(result.final_state.biomass[0].shape)[0],
        weights=result.final_state.biomass[0],
    )
    assert final_center < initial_center


def test_chemotaxis_requires_a_declared_field():
    with pytest.raises(EcosystemError, match="unknown field"):
        config(
            species=(
                SpeciesConfig(
                    "invalid",
                    0.1,
                    0.1,
                    (0.2, 0.2),
                    (1.0, 1.0),
                    chemotaxis_field="missing",
                    chemotaxis_sensitivity=1.0,
                ),
                config().species[1],
            )
        )


def test_quorum_switching_activates_a_state_with_hysteresis_and_dwell():
    result = run(
        config(
            duration_h=0.2,
            timestep_h=0.01,
            nutrients=(NutrientConfig("food", 1.0, 0.0),),
            species=(
                SpeciesConfig(
                    "switcher",
                    0.0,
                    0.5,
                    (0.2,),
                    (1.0,),
                    seed_regions=(SeedRegion(0.5, 0.5, 0.2, 0.8),),
                    phenotypes=(
                        PhenotypeConfig(
                            "biofilm",
                            activation_threshold=0.1,
                            deactivation_threshold=0.05,
                            growth_multiplier=0.2,
                            minimum_dwell_h=0.05,
                        ),
                    ),
                ),
            ),
        )
    )
    assert np.any(result.final_state.phenotype_indices[0] == 1)
    assert result.frames[-1].phenotype_indices is not None
    assert result.frames[-1].to_dict(("switcher",), ("food",), 1.0)["phenotypes"]


def test_phenotype_thresholds_are_validated():
    with pytest.raises(EcosystemError, match="deactivation_threshold"):
        PhenotypeConfig("invalid", 0.2, 0.4)


def test_spatial_capability_rate_and_limiting_factor_are_exported():
    result = run(
        config(
            duration_h=0.1,
            nutrients=(NutrientConfig("carbon", 1.0, 0.0),),
            conditions=(
                ConditionConfig(
                    "oxygen",
                    0.0,
                    10.0,
                    boundary_value=1.0,
                    boundary_edges=("top",),
                ),
            ),
            species=(
                SpeciesConfig(
                    "aerobe",
                    0.1,
                    1.0,
                    (0.2,),
                    (1.0,),
                    capabilities=(
                        Capability(
                            "aerobic-respiration",
                            1.0,
                            "carbon",
                            0.2,
                            oxygen_half_saturation=0.1,
                        ),
                    ),
                ),
            ),
        )
    )
    frame = result.frames[-1]
    assert frame.niche_rates[0][0, 0] > frame.niche_rates[0][-1, 0]
    assert frame.niche_limiting_factors[0][-1, 0] == "oxygen"
    exported = frame.to_dict(("aerobe",), ("carbon",), 1.0, ("oxygen",))
    assert exported["niche"]["limiting_factor"]["aerobe"][-1][0] == "oxygen"


def test_hill_additive_response_supports_inhibition_and_stimulation():
    concentrations = np.array([0.0, 1.0, 10.0])
    inhibition = apply_effect(
        AdditiveEffect("drug", 0.1, 1.0, 1.0, direction="decreasing"), concentrations
    )
    stimulation = apply_effect(AdditiveEffect("signal", 1.0, 2.0, 1.0), concentrations)
    assert hill_response(concentrations, 1.0, 1.0)[0] == 0.0
    assert inhibition[0] == pytest.approx(1.0)
    assert inhibition[-1] < inhibition[0]
    assert stimulation[-1] > stimulation[0]


def test_additive_field_diffuses_decays_and_changes_growth(tmp_path: Path):
    result = run(
        config(
            duration_h=0.2,
            nutrients=(NutrientConfig("carbon", 1.0, 0.0),),
            additives=(
                AdditiveConfig(
                    "drug",
                    1.0,
                    10.0,
                    decay_per_h=1.0,
                    boundary_value=0.0,
                    boundary_edges=("bottom",),
                ),
            ),
            species=(
                SpeciesConfig(
                    "target",
                    0.1,
                    1.0,
                    (0.2,),
                    (1.0,),
                    additive_effects=(AdditiveEffect("drug", 0.1, 1.0, 0.2, 2.0),),
                ),
            ),
        )
    )
    assert result.final_state.additives[0].mean() < result.frames[0].additives[0].mean()
    assert (
        result.final_state.biomass[0].mean()
        < run(
            config(
                duration_h=0.2,
                nutrients=(NutrientConfig("carbon", 1.0, 0.0),),
                species=(SpeciesConfig("target", 0.1, 1.0, (0.2,), (1.0,)),),
            )
        )
        .final_state.biomass[0]
        .mean()
    )
    exported = result.write_frames(tmp_path / "frames.json")
    assert '"additives"' in exported.read_text()
