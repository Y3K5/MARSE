"""Oxygen roles: how a species lives with oxygen, checked against its processes.

Stage 2d, increment 2d.1 (docs/stage-2d-plan.md). Every species declares one
of five roles, and a network whose processes break its species' roles is
refused as it loads (criterion D1). An obligate anaerobe declared correctly
grows without oxygen at the rate its other factors give (D2): the requirement
behind the version 1 engine's known defect 3, on the version 2 engine. The
examples gain roles without any change to their results (D3); the digests
before and after are recorded in docs/validation.md.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from marse.cli import main
from marse.core.config import ConfigError
from marse.core.well_mixed import run
from marse.schemas import OXYGEN_ROLES, experiment_from_dict, load_network, network_from_dict

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = [
    ROOT / "examples" / "networks" / "glucose_cross_feeding.json",
    ROOT / "examples" / "networks" / "surface_biofilm_1d.json",
    ROOT / "examples" / "networks" / "surface_biofilm_3d.json",
    ROOT / "examples" / "environments" / "dental" / "dental_surfaces.json",
]

AEROBIC = ["oxygen", "carbon_dioxide", "ammonium"]
LACTIC = ["lactate", "carbon_dioxide", "ammonium"]


def monod(component: str, k: float = 0.05) -> dict:
    return {"component": component, "form": "monod", "half_saturation_mol_per_m3": k}


def inhibition(component: str, k: float = 0.01) -> dict:
    return {"component": component, "form": "inhibition", "inhibition_mol_per_m3": k}


def haldane(component: str) -> dict:
    return {
        "component": component,
        "form": "haldane",
        "half_saturation_mol_per_m3": 0.005,
        "inhibition_mol_per_m3": 0.05,
    }


def growth(name: str, balanced_by: list[str], *oxygen_factors: dict) -> dict:
    return {
        "name": name,
        "kind": "growth",
        "biomass": "bug",
        "substrate": "glucose",
        "yield_mol_per_mol": 3.0 if balanced_by is AEROBIC else 0.6,
        "balanced_by": balanced_by,
        "rate": {
            "maximum_per_h": 0.4,
            "factors": [monod("glucose"), monod("ammonium"), *oxygen_factors],
        },
    }


def respire(*oxygen: dict) -> dict:
    """Growth that consumes oxygen; it needs a factor on oxygen, a monod one by default."""
    return growth("respire", AEROBIC, *(oxygen or (monod("oxygen", 0.005),)))


def ferment(*oxygen: dict) -> dict:
    """Growth that makes lactate and consumes no oxygen."""
    return growth("ferment", LACTIC, *oxygen)


ENDOGENOUS = {
    "name": "endogenous",
    "kind": "reaction",
    "stoichiometry_mol_per_mol": {"bug": -1},
    "balanced_by": AEROBIC,
    "rate": {"maximum_per_h": 0.01, "proportional_to": "bug", "factors": [monod("oxygen")]},
}


def raw_network(role: str | None, *processes: dict, oxygen: bool = True) -> dict:
    components = [
        {"name": "glucose", "phase": "dissolved", "formula": "C6H12O6"},
        {"name": "ammonium", "phase": "dissolved", "formula": "NH4", "charge": 1},
        {"name": "carbon_dioxide", "phase": "dissolved", "formula": "CO2"},
        {"name": "lactate", "phase": "dissolved", "formula": "C3H5O3", "charge": -1},
        {"name": "bug", "phase": "particulate", "formula": "CH1.8O0.5N0.2"},
    ]
    if oxygen:
        components.insert(1, {"name": "oxygen", "phase": "dissolved", "formula": "O2"})
    if role is not None:
        components[-1]["oxygen_role"] = role
    return {"schema_version": 2, "components": components, "processes": list(processes)}


# --- D1: what each role accepts, and what it refuses ------------------------------------


@pytest.mark.parametrize(
    ("role", "processes"),
    [
        ("obligate_aerobe", [respire()]),
        ("obligate_aerobe", [respire(), ENDOGENOUS]),
        ("microaerophile", [respire(haldane("oxygen"))]),
        ("facultative", [respire(), ferment()]),
        ("facultative", [respire(), ferment(inhibition("oxygen"))]),
        ("aerotolerant", [ferment()]),
        ("aerotolerant", [ferment(inhibition("oxygen"))]),
        ("obligate_anaerobe", [ferment(inhibition("oxygen"))]),
    ],
)
def test_each_role_accepts_processes_that_obey_it(role, processes):
    loaded = network_from_dict(raw_network(role, *processes))
    assert loaded.component("bug").oxygen_role == role


def test_an_obligate_anaerobe_in_a_network_without_oxygen_needs_no_inhibition():
    loaded = network_from_dict(raw_network("obligate_anaerobe", ferment(), oxygen=False))
    assert loaded.component("bug").oxygen_role == "obligate_anaerobe"


@pytest.mark.parametrize(
    ("role", "processes", "message"),
    [
        (
            "obligate_aerobe",
            [respire(), ferment()],
            "is an obligate aerobe, but its growth process 'ferment' does not consume oxygen",
        ),
        (
            "microaerophile",
            [ferment()],
            "is a microaerophile, but its growth process 'ferment' does not consume oxygen",
        ),
        (
            "microaerophile",
            [respire()],
            "its growth process 'respire' has no haldane factor on oxygen",
        ),
        (
            "facultative",
            [respire()],
            "is a facultative species, but every one of its growth processes consumes oxygen",
        ),
        (
            "facultative",
            [ferment()],
            "is a facultative species, but none of its growth processes consumes oxygen",
        ),
        (
            "aerotolerant",
            [respire()],
            "is an aerotolerant species, but its process 'respire' consumes oxygen",
        ),
        (
            "aerotolerant",
            [ferment(), ENDOGENOUS],
            "is an aerotolerant species, but its process 'endogenous' consumes oxygen",
        ),
        (
            "aerotolerant",
            [ferment(monod("oxygen"))],
            "its process 'ferment' has a monod factor on 'oxygen', so it could not run "
            "without oxygen",
        ),
        (
            "obligate_anaerobe",
            [ferment(haldane("oxygen"))],
            "its process 'ferment' has a haldane factor on 'oxygen'",
        ),
        (
            "obligate_anaerobe",
            [respire()],
            "is an obligate anaerobe, but its process 'respire' consumes oxygen",
        ),
        (
            "obligate_anaerobe",
            [ferment()],
            "its growth process 'ferment' has no inhibition factor on oxygen",
        ),
    ],
)
def test_each_role_refuses_processes_that_break_it_naming_the_process(role, processes, message):
    with pytest.raises(ConfigError, match=message) as refused:
        network_from_dict(raw_network(role, *processes))
    assert "network.components 'bug'" in str(refused.value)


@pytest.mark.parametrize("role", ["obligate_aerobe", "microaerophile", "facultative"])
def test_a_species_that_grows_on_oxygen_needs_oxygen_in_the_network(role):
    with pytest.raises(ConfigError, match="no component of this network is oxygen"):
        network_from_dict(raw_network(role, ferment(), oxygen=False))


def test_the_version_one_periodontal_error_cannot_be_written_down():
    # Known defect 3: anaerobes whose growth had a Monod term on oxygen, so
    # that without oxygen they could not grow at all.
    raw = raw_network("obligate_anaerobe", ferment(monod("oxygen"), inhibition("oxygen")))
    with pytest.raises(ConfigError, match="appears more than once"):
        network_from_dict(raw)  # one factor per component already refuses both together
    with pytest.raises(ConfigError, match="has a monod factor on 'oxygen'"):
        network_from_dict(raw_network("obligate_anaerobe", ferment(monod("oxygen"))))


def test_a_species_without_a_role_is_refused_naming_its_growth_process():
    with pytest.raises(
        ConfigError, match="'bug': grows in process 'respire', so it is a species and needs"
    ):
        network_from_dict(raw_network(None, respire()))


def test_a_particulate_component_that_does_not_grow_needs_no_role():
    loaded = network_from_dict(raw_network(None, oxygen=True))
    assert loaded.component("bug").oxygen_role is None


def test_a_role_on_a_dissolved_component_is_refused():
    raw = raw_network("obligate_aerobe", respire())
    raw["components"][0]["oxygen_role"] = "obligate_aerobe"
    with pytest.raises(ConfigError, match="'glucose' is dissolved; an oxygen role belongs to"):
        network_from_dict(raw)


def test_an_unknown_role_is_refused_listing_the_roles():
    with pytest.raises(ConfigError, match="expected one of 'obligate_aerobe'"):
        network_from_dict(raw_network("anaerobe", ferment()))


def test_oxygen_is_found_by_its_formula_not_its_name():
    raw = raw_network("obligate_aerobe", respire())
    text = json.dumps(raw).replace('"oxygen"', '"o2"')
    loaded = network_from_dict(json.loads(text))
    assert loaded.component("o2").is_oxygen
    # Ozone is not oxygen, and neither is peroxide, a charged O2.
    for formula, charge in (("O3", 0), ("O2", -2)):
        raw = raw_network("aerotolerant", ferment())
        raw["components"][1].update(formula=formula, charge=charge)
        assert not network_from_dict(raw).component("oxygen").is_oxygen


def test_the_role_survives_the_trip_through_json():
    loaded = network_from_dict(raw_network("facultative", respire(), ferment()))
    written = loaded.to_dict()
    assert written["components"][-1]["oxygen_role"] == "facultative"
    assert network_from_dict(json.loads(json.dumps(written))) == loaded


def test_every_role_is_documented():
    text = (ROOT / "docs" / "networks.md").read_text("utf-8")
    for role in OXYGEN_ROLES:
        assert f"| `{role}` |" in text


# --- D2: an obligate anaerobe grows without oxygen -------------------------------------


def anoxic_culture(with_oxygen: bool) -> dict:
    """An obligate anaerobe fermenting in a closed box that holds no oxygen."""
    raw = raw_network("obligate_anaerobe", ferment(*([inhibition("oxygen")] * with_oxygen)))
    if not with_oxygen:
        raw["components"] = [c for c in raw["components"] if c["name"] != "oxygen"]
    return raw | {
        "experiment_id": "anoxic",
        "initial_mol_per_m3": {"glucose": 5.0, "ammonium": 2.0, "bug": 0.05},
        "duration_h": 12.0,
        "timestep_h": 0.5,
        "relative_tolerance": 1e-9,
    }


def test_an_obligate_anaerobe_grows_without_oxygen_as_its_other_factors_say():
    # With no oxygen, its inhibition factor is exactly 1, so the culture must
    # grow exactly as the same culture in a network without oxygen at all: bit
    # for bit, since multiplying by 1 is exact and oxygen's error is zero.
    inhibited = run(experiment_from_dict(anoxic_culture(with_oxygen=True)))
    alone = run(experiment_from_dict(anoxic_culture(with_oxygen=False)))
    grown = inhibited.final_mol_per_m3["bug"] / 0.05
    assert grown > 10, f"an anaerobe without oxygen grew only {grown:.2f}-fold"
    columns = [inhibited.component_names.index(n) for n in alone.component_names]
    np.testing.assert_array_equal(
        inhibited.concentrations_mol_per_m3[:, columns], alone.concentrations_mol_per_m3
    )
    oxygen = inhibited.component_names.index("oxygen")
    assert np.all(inhibited.concentrations_mol_per_m3[:, oxygen] == 0)


def test_oxygen_slows_an_obligate_anaerobe_through_its_inhibition_factor():
    # Oxygen held at the inhibition constant halves the anaerobe's rate; at
    # forty times it, a forty-first is left. Nothing consumes the oxygen, and
    # the step is short enough for the growth over it to measure the rate.
    def grown(oxygen: float) -> float:
        raw = anoxic_culture(with_oxygen=True)
        raw["initial_mol_per_m3"]["oxygen"] = oxygen
        raw["duration_h"], raw["timestep_h"] = 1e-4, 1e-4
        result = run(experiment_from_dict(raw))
        return result.final_mol_per_m3["bug"] - 0.05

    anoxic = grown(0.0)
    assert grown(0.01) / anoxic == pytest.approx(1 / 2, rel=1e-4)
    assert grown(0.4) / anoxic == pytest.approx(1 / 41, rel=1e-4)


# --- D3: the examples gain roles; marse check reports them ------------------------------


@pytest.mark.parametrize("path", EXAMPLES, ids=lambda p: p.name)
def test_every_example_species_has_a_role(path):
    network = load_network(path) if "networks" in path.parts else None
    if network is None:
        network = experiment_from_dict(json.loads(path.read_text("utf-8"))).network
    grown = {p.growth.biomass for p in network.processes if p.growth is not None}
    assert grown
    assert all(network.component(name).oxygen_role in OXYGEN_ROLES for name in grown)


def test_marse_check_prints_each_species_role(capsys):
    assert main(["check", str(EXAMPLES[0])]) == 0
    out = capsys.readouterr().out
    assert "oxygen roles, each checked against its processes: heterotroph obligate_aerobe," in out
    assert "fermenter aerotolerant" in out
