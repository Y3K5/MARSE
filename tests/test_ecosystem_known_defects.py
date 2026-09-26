"""Known defects in the ecosystem engine, written as the behaviour it must have.

Every test here states a physical or reproducibility requirement that the
ecosystem engine currently breaks. Each defect was found in a review of the
engine and confirmed by running it; ``docs/validation.md`` lists them with the
evidence. They are marked ``xfail(strict=True, raises=AssertionError)``:

- while the defect exists, the requirement's assertion fails and the test is
  reported as an expected failure, so the suite stays honest about it;
- the fix is what flips it. A passing test is then reported as a failure, so
  the marker has to be removed in the same change and the requirement becomes
  an ordinary regression test;
- any error other than the requirement's own assertion (a crash, a refused
  configuration) fails the test outright, so an xfail can never hide one.

These are new tests that describe unfinished work. None of them replaces or
disables an existing test.
"""

import json
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from marse.cli import main
from marse.ecosystem import (
    ConditionConfig,
    EcosystemConfig,
    EcosystemError,
    NutrientConfig,
    SeedRegion,
    SpeciesConfig,
    run,
)
from marse.ecosystem.model import ecosystem_from_dict
from marse.microbes.niche import Capability
from marse.spatial.solutes import oxygen_diffusivity_um2_per_h

EXPERIMENTS = Path(__file__).resolve().parents[1] / "examples" / "experiments"


def known_defect(number: int, summary: str, fixed_in: str):
    return pytest.mark.xfail(
        strict=True,
        raises=AssertionError,
        reason=f"known defect {number} (docs/validation.md): {summary}; fixed in {fixed_in}",
    )


def well_mixed(**overrides) -> EcosystemConfig:
    """A 2x2 closed box with uniform fields, so transport has nothing to do."""
    values = {
        "experiment_id": "known-defect",
        "width": 2,
        "height": 2,
        "cell_size_um": 10.0,
        "duration_h": 1.0,
        "timestep_h": 0.01,
        "seed": 1,
        "carrying_capacity": 1e6,
        "nutrients": (NutrientConfig("carbon", 2.0, 10.0),),
        "species": (SpeciesConfig("grower", 0.01, 0.5, (0.2,), (1.0,)),),
    }
    return EcosystemConfig(**(values | overrides))


def totals(result, frame: int = -1) -> tuple[np.ndarray, np.ndarray]:
    state = result.frames[frame]
    return state.biomass.sum(axis=(1, 2)), state.nutrients.sum(axis=(1, 2))


@known_defect(2, "uptake follows potential growth, not actual growth", "Stage 2")
def test_no_substrate_is_consumed_when_growth_is_zero():
    # Biomass starts at the carrying capacity, so the logistic term forbids any
    # growth. With no maintenance declared, nothing may be consumed either.
    config = well_mixed(
        carrying_capacity=1.0,
        species=(SpeciesConfig("grower", 1.0, 0.5, (0.2,), (1.0,)),),
    )
    result = run(config)
    biomass_before, carbon_before = totals(result, 0)
    biomass_after, carbon_after = totals(result)
    np.testing.assert_allclose(biomass_after, biomass_before, rtol=1e-12)
    np.testing.assert_allclose(carbon_after, carbon_before, rtol=1e-12)


@known_defect(2, "consumption is not growth divided by yield", "Stage 2")
def test_each_nutrient_is_consumed_at_growth_over_its_yield():
    # Two required nutrients with different yields. The second is scarce
    # relative to its half-saturation, so it limits growth strongly.
    config = well_mixed(
        nutrients=(NutrientConfig("carbon", 1.0, 10.0), NutrientConfig("nitrogen", 2.0, 10.0)),
        species=(SpeciesConfig("grower", 0.01, 0.5, (0.1, 20.0), (0.5, 2.0)),),
    )
    result = run(config)
    biomass_before, nutrients_before = totals(result, 0)
    biomass_after, nutrients_after = totals(result)
    grown = float(biomass_after.sum() - biomass_before.sum())
    consumed = nutrients_before - nutrients_after
    np.testing.assert_allclose(consumed, grown / np.array([0.5, 2.0]), rtol=1e-9)


@known_defect(2, "production adds material that no consumed substrate pays for", "Stage 2")
def test_production_does_not_create_matter():
    # One unit of carbon makes one unit of biomass (yield 1), and the species
    # also releases a product. The product can only come out of the carbon, so
    # the total of carbon, product and biomass must not rise.
    config = well_mixed(
        nutrients=(NutrientConfig("carbon", 2.0, 10.0), NutrientConfig("product", 1.0, 10.0)),
        species=(
            SpeciesConfig(
                "grower",
                0.01,
                0.5,
                (0.2, 1e-9),  # the product is never limiting
                (1.0, 1e12),  # and is not consumed in any measurable amount
                production_per_nutrient=(0.0, 0.5),
            ),
        ),
    )
    result = run(config)
    before = sum(float(x.sum()) for x in totals(result, 0))
    after = sum(float(x.sum()) for x in totals(result))
    assert after <= before * (1.0 + 1e-12), f"total rose from {before!r} to {after!r}"


@known_defect(4, "a capability repeats the substrate Monod term", "Stage 2")
def test_a_capability_does_not_apply_the_substrate_limit_twice():
    half_saturation = 0.35
    capability = Capability("carbon-use", 0.5, "carbon", half_saturation)
    config = well_mixed(
        duration_h=1e-3,
        timestep_h=1e-3,
        nutrients=(NutrientConfig("carbon", half_saturation, 10.0),),
        species=(
            SpeciesConfig(
                "grower", 0.01, 0.5, (half_saturation,), (1.0,), capabilities=(capability,)
            ),
        ),
    )
    result = run(config)
    biomass_before, _ = totals(result, 0)
    biomass_after, _ = totals(result)
    rate = math.log(biomass_after[0] / biomass_before[0]) / config.timestep_h
    # At C = K a single Monod term halves the maximum rate; applied twice it quarters it.
    assert rate == pytest.approx(0.5 * 0.5, rel=1e-3)


@known_defect(7, "species are updated one after another within a step", "Stage 2")
def test_the_outcome_does_not_depend_on_the_order_species_are_listed():
    first = SpeciesConfig("fast", 0.2, 0.9, (0.5,), (1.0,))
    second = SpeciesConfig("thrifty", 0.2, 0.5, (0.05,), (1.0,))
    config = well_mixed(
        duration_h=2.0,
        timestep_h=0.05,
        nutrients=(NutrientConfig("carbon", 0.3, 10.0),),
        species=(first, second),
    )
    listed = run(config).final_state.biomass.sum(axis=(1, 2))
    swapped = run(replace(config, species=(second, first))).final_state.biomass.sum(axis=(1, 2))
    np.testing.assert_allclose(swapped[::-1], listed, rtol=1e-12)


@known_defect(8, "each species has its own carrying capacity", "Stage 2")
def test_species_share_one_carrying_capacity():
    config = well_mixed(
        duration_h=40.0,
        timestep_h=0.05,
        carrying_capacity=1.0,
        nutrients=(NutrientConfig("carbon", 1e6, 10.0),),
        species=(
            SpeciesConfig("one", 0.1, 0.8, (0.01,), (1e6,)),
            SpeciesConfig("two", 0.1, 0.8, (0.01,), (1e6,)),
        ),
    )
    per_cell = run(config).final_state.biomass.sum(axis=0)
    assert per_cell.max() <= config.carrying_capacity * (1.0 + 1e-9), (
        f"{per_cell.max():.3f} per cell against a capacity of {config.carrying_capacity}"
    )


@known_defect(5, "negative values are clipped, which creates biomass", "Stage 2")
def test_chemotaxis_conserves_biomass_or_refuses_the_step():
    # A fixed, steep attractant gradient on the top edge and no growth at all:
    # chemotaxis only moves biomass, so the total must not change. A timestep
    # too large for the advection has to be refused, not clipped.
    config = well_mixed(
        width=4,
        height=4,
        duration_h=0.5,
        timestep_h=0.5,
        conditions=(ConditionConfig("attractant", 0.0, 0.0, 100.0, ("top",)),),
        species=(
            SpeciesConfig(
                "swimmer",
                1.0,
                0.0,
                (0.2,),
                (1.0,),
                chemotaxis_field="attractant",
                chemotaxis_sensitivity=10.0,
            ),
        ),
    )
    try:
        result = run(config)
    except EcosystemError:
        return  # refused: the correct outcome for an unstable step
    before, _ = totals(result, 0)
    after, _ = totals(result)
    np.testing.assert_allclose(after, before, rtol=1e-12)


@known_defect(6, "mutation is a property of grid cells, not of lineages", "Stage 2")
def test_empty_space_cannot_carry_a_mutation():
    colony = SeedRegion(0.0, 0.0, 0.3, 0.5)
    config = well_mixed(
        width=6,
        height=6,
        mutation_interval_h=0.01,
        duration_h=0.02,
        species=(
            SpeciesConfig(
                "grower",
                0.0,
                0.5,
                (0.2,),
                (1.0,),
                mutation_probability=1.0,
                seed_regions=(colony,),
            ),
        ),
    )
    result = run(config)
    empty = result.frames[0].biomass[0] == 0.0
    assert empty.any(), "the test needs empty cells to be meaningful"
    assert not result.final_state.mutations[0][empty].any(), "empty cells were marked mutant"


@known_defect(3, "the periodontal anaerobes need oxygen to grow", "Stage 2 (O2 role)")
def test_the_periodontal_anaerobes_grow_without_oxygen():
    raw = json.loads((EXPERIMENTS / "periodontal_pathogen_biofilm.json").read_text("utf-8"))
    raw["duration_h"] = 2.0
    for field in raw["nutrients"]:
        if field["name"] == "oxygen":
            field["initial"] = 0.0
            field["boundary_value"] = None
            field["boundary_edges"] = []
    result = run(ecosystem_from_dict(raw))
    before, _ = totals(result, 0)
    after, _ = totals(result)
    fold = after / before
    names = [species["name"] for species in raw["species"]]
    assert np.all(fold > 1.01), f"growth without oxygen: {dict(zip(names, fold, strict=True))}"


@known_defect(1, "example oxygen diffusivities are about 4e5 times too small", "Stage 2")
@pytest.mark.parametrize(
    "experiment", ["two_species_ecosystem.json", "periodontal_pathogen_biofilm.json"]
)
def test_example_oxygen_diffusivity_is_physical(experiment: str):
    raw = json.loads((EXPERIMENTS / experiment).read_text("utf-8"))
    temperature = next(
        (c["initial"] for c in raw.get("conditions", []) if c["name"] == "temperature_c"), 25.0
    )
    oxygen = next(n for n in raw["nutrients"] if n["name"] == "oxygen")
    in_water = oxygen_diffusivity_um2_per_h(temperature)
    # Effective diffusivity in a biofilm lies below that in water, but not by
    # more than an order of magnitude for a small solute such as oxygen.
    assert 0.1 * in_water <= oxygen["diffusivity"] <= in_water, (
        f"{oxygen['diffusivity']} um2/h against {in_water:.3g} um2/h in water"
    )


# Fixed in Stage 1: ecosystem runs write a manifest and replay exactly. This is
# now an ordinary regression test (tests/test_ecosystem_provenance.py goes further).
def test_an_ecosystem_run_writes_a_replayable_manifest(tmp_path: Path):
    experiment = {
        "experiment_id": "manifest-check",
        "width": 2,
        "height": 2,
        "cell_size_um": 10.0,
        "duration_h": 0.05,
        "timestep_h": 0.01,
        "seed": 3,
        "nutrients": [{"name": "carbon", "initial": 1.0, "diffusivity": 10.0}],
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
    path = tmp_path / "experiment.json"
    path.write_text(json.dumps(experiment), encoding="utf-8")
    output = tmp_path / "run"
    assert main(["ecosystem", str(path), "-o", str(output)]) == 0
    manifest = output / "manifest.json"
    assert manifest.is_file(), "no manifest.json next to the ecosystem outputs"
    assert main(["replay", str(manifest)]) == 0
