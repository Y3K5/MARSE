"""Reaction networks: every process is proven to conserve carbon, nitrogen and electrons.

The rows MARSE derives are checked against textbook reactions and, for growth,
against the half-reaction method of Rittmann and McCarty (2001, chapter 2),
an independent route to the same stoichiometry. The remaining tests pin down
what the loader refuses and what its messages say, because a refused process
is useful only if the message says what to change.
"""

import json
import random
import re
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import pytest

from marse.cli import main
from marse.core.config import ConfigError
from marse.schemas import ContinuityError, load_network, network_from_dict
from marse.schemas._reading import DIMENSIONLESS, NUMERIC_KINDS, UNIT_SUFFIXES
from marse.schemas.network import SCHEMA

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "networks" / "glucose_cross_feeding.json"

COMPONENTS = {
    "glucose": ("dissolved", "C6H12O6", 0),
    "oxygen": ("dissolved", "O2", 0),
    "ammonium": ("dissolved", "NH4", 1),
    "carbon_dioxide": ("dissolved", "CO2", 0),
    "nitrate": ("dissolved", "NO3", -1),
    "ethanol": ("dissolved", "C2H6O", 0),
    "lactate": ("dissolved", "C3H5O3", -1),
    "heterotroph": ("particulate", "CH1.8O0.5N0.2", 0),
    "cells": ("particulate", "C5H7O2N", 0),
}
CLOSE_AEROBICALLY = ["oxygen", "carbon_dioxide", "ammonium"]


def component(name: str) -> dict:
    phase, formula, charge = COMPONENTS[name]
    return {"name": name, "phase": phase, "formula": formula, "charge": charge}


def raw_network(*processes: dict) -> dict:
    return {
        "schema_version": 2,
        "components": [component(name) for name in COMPONENTS],
        "processes": list(processes),
    }


def network(*processes: dict):
    return network_from_dict(raw_network(*processes))


def reaction(name: str = "reaction", balanced_by=(), **stoichiometry) -> dict:
    return {
        "name": name,
        "kind": "reaction",
        "stoichiometry_mol_per_mol": stoichiometry,
        "balanced_by": list(balanced_by),
    }


def growth(biomass, substrate, yield_, balanced_by, name="growth", products=None) -> dict:
    raw = {
        "name": name,
        "kind": "growth",
        "biomass": biomass,
        "substrate": substrate,
        "yield_mol_per_mol": yield_,
        "balanced_by": list(balanced_by),
    }
    if products is not None:
        raw["products_mol_per_mol"] = products
    return raw


def only_process(*processes: dict):
    return network(*processes).processes[0]


# --- textbook chemistry ------------------------------------------------------


def test_respiration_of_glucose():
    # C6H12O6 + 6 O2 -> 6 CO2 + 6 H2O
    p = only_process(reaction(glucose=-1, balanced_by=["oxygen", "carbon_dioxide"]))
    assert p.coefficients == {"glucose": -1, "oxygen": -6, "carbon_dioxide": 6}
    assert (p.water, p.protons) == (6, 0)


def test_alcoholic_fermentation():
    # Gay-Lussac: C6H12O6 -> 2 C2H5OH + 2 CO2
    p = only_process(reaction(glucose=-1, balanced_by=["ethanol", "carbon_dioxide"]))
    assert p.coefficients == {"glucose": -1, "carbon_dioxide": 2, "ethanol": 2}
    assert (p.water, p.protons) == (0, 0)


def test_homolactic_fermentation_balances_carbon_dioxide_to_zero():
    # C6H12O6 -> 2 CH3CHOHCOO- + 2 H+
    p = only_process(reaction(glucose=-1, balanced_by=["lactate", "carbon_dioxide"]))
    assert p.coefficients == {"glucose": -1, "lactate": 2}
    assert (p.water, p.protons) == (0, 2)


def test_nitrification():
    # NH4+ + 2 O2 -> NO3- + H2O + 2 H+
    p = only_process(reaction(ammonium=-1, balanced_by=["oxygen", "nitrate"]))
    assert p.coefficients == {"oxygen": -2, "ammonium": -1, "nitrate": 1}
    assert (p.water, p.protons) == (1, 2)


def test_the_worked_example_of_aerobic_growth():
    # docs/networks.md works this one through by hand.
    p = only_process(growth("heterotroph", "glucose", 3.6, CLOSE_AEROBICALLY))
    assert p.coefficients == {
        "glucose": F(-5, 18),
        "oxygen": F(-37, 60),
        "ammonium": F(-1, 5),
        "carbon_dioxide": F(2, 3),
        "heterotroph": 1,
    }
    assert (p.water, p.protons) == (F(16, 15), F(1, 5))
    assert p.equation() == (
        "0.277778 glucose + 0.616667 oxygen + 0.2 ammonium -> "
        "0.666667 carbon_dioxide + heterotroph + 1.06667 H2O + 0.2 H+"
    )


def test_endogenous_respiration_uses_a_quarter_mol_of_oxygen_per_electron():
    p = only_process(reaction(heterotroph=-1, balanced_by=CLOSE_AEROBICALLY))
    assert p.coefficient("oxygen") == F(-21, 20)  # 4.2 electrons / 4
    assert p.coefficient("ammonium") == F(1, 5)  # the nitrogen comes back out


# --- the half-reaction method (Rittmann and McCarty 2001, chapter 2) ---------

# Half reactions per electron equivalent, all written as reductions; negative
# coefficients are on the left. Rc supplies part of the cell carbon as bicarbonate.
RA_OXYGEN = {"O2": F(-1, 4), "H+": -1, "e-": -1, "H2O": F(1, 2)}
RC_CELLS_ON_AMMONIUM = {
    "CO2": F(-1, 5),
    "HCO3-": F(-1, 20),
    "NH4+": F(-1, 20),
    "H+": -1,
    "e-": -1,
    "C5H7O2N": F(1, 20),
    "H2O": F(9, 20),
}
RD_GLUCOSE = {"CO2": F(-1, 4), "H+": -1, "e-": -1, "C6H12O6": F(1, 24), "H2O": F(1, 4)}

ATOMS = {
    "O2": ({"O": 2}, 0),
    "H+": ({"H": 1}, 1),
    "e-": ({}, -1),
    "H2O": ({"H": 2, "O": 1}, 0),
    "CO2": ({"C": 1, "O": 2}, 0),
    "HCO3-": ({"H": 1, "C": 1, "O": 3}, -1),
    "NH4+": ({"N": 1, "H": 4}, 1),
    "C5H7O2N": ({"C": 5, "H": 7, "O": 2, "N": 1}, 0),
    "C6H12O6": ({"C": 6, "H": 12, "O": 6}, 0),
}


@pytest.mark.parametrize("half", [RA_OXYGEN, RC_CELLS_ON_AMMONIUM, RD_GLUCOSE])
def test_the_transcribed_half_reactions_balance(half):
    for element in "CHON":
        assert sum(c * ATOMS[s][0].get(element, 0) for s, c in half.items()) == 0, element
    assert sum(c * ATOMS[s][1] for s, c in half.items()) == 0, "charge"


@pytest.mark.parametrize("fs", [F(1, 2), F(3, 5), F(7, 10)])
def test_growth_rows_reproduce_the_half_reaction_method(fs):
    # R = fe Ra + fs Rc - Rd, with fe + fs = 1.
    overall: dict[str, F] = {}
    for half, weight in ((RA_OXYGEN, 1 - fs), (RC_CELLS_ON_AMMONIUM, fs), (RD_GLUCOSE, F(-1))):
        for species, c in half.items():
            overall[species] = overall.get(species, F(0)) + weight * c
    assert overall.pop("e-") == 0
    bicarbonate = overall.pop("HCO3-")  # HCO3- + H+ -> CO2 + H2O
    overall["CO2"] += bicarbonate
    overall["H2O"] += bicarbonate
    overall["H+"] -= bicarbonate
    per_cell = {species: c / overall["C5H7O2N"] for species, c in overall.items()}

    yield_ = 1 / -per_cell["C6H12O6"]  # mol of cells per mol of glucose
    assert F(repr(float(yield_))) == yield_  # a decimal, so the file states it exactly
    p = only_process(growth("cells", "glucose", float(yield_), CLOSE_AEROBICALLY))
    assert p.coefficients == {
        "glucose": per_cell["C6H12O6"],
        "oxygen": per_cell["O2"],
        "ammonium": per_cell["NH4+"],
        "carbon_dioxide": per_cell["CO2"],
        "cells": 1,
    }
    assert (p.water, p.protons) == (per_cell["H2O"], per_cell["H+"])


# --- what a process may not do ------------------------------------------------


def test_production_from_nothing_is_refused_naming_the_process_and_quantities():
    # Known defect 2's pattern: biomass and a product appear that the one mol of
    # glucose consumed cannot pay for.
    with pytest.raises(ContinuityError) as refused:
        network(reaction("free_lunch", glucose=-1, heterotroph=6, lactate=1))
    assert refused.value.process == "free_lunch"
    assert refused.value.imbalance == {"carbon": 3, "nitrogen": F(6, 5), "electrons": F(66, 5)}
    assert "'free_lunch': creates 3 mol C, 1.2 mol N and 13.2 mol e- per unit" in str(refused.value)


def test_destruction_is_refused_too():
    with pytest.raises(
        ContinuityError, match=re.escape("destroys 1 mol C, 0.2 mol N and 4.2 mol e-")
    ):
        network(reaction("vanishing", heterotroph=-1))


def test_a_row_that_nearly_balances_is_still_refused():
    # A coefficient rounded when it was copied: exact balance means exact.
    with pytest.raises(ContinuityError, match=r"destroys 0\.001 mol C per unit of reaction"):
        network(reaction(glucose=-1, oxygen=-6, carbon_dioxide=5.999))


def test_a_product_the_electrons_cannot_pay_for_is_refused():
    with pytest.raises(
        ContinuityError, match=re.escape("0.2 mol e- per mol of heterotroph formed")
    ):
        network(
            growth(
                "heterotroph",
                "glucose",
                0.6,
                ["carbon_dioxide", "ammonium"],
                products={"lactate": 1.8},
            )
        )


def test_a_quantity_nothing_in_balanced_by_carries_is_named():
    with pytest.raises(
        ContinuityError,
        match="nothing in balanced_by carries nitrogen; add a component that does: "
        "the nitrogen source",
    ):
        network(growth("heterotroph", "glucose", 3.6, ["oxygen", "carbon_dioxide"]))


def test_closers_that_cannot_balance_everything_together_are_refused():
    # One biomass cannot absorb glucose's carbon and electrons in the right ratio.
    with pytest.raises(ContinuityError, match="no choice of coefficients for 'heterotroph'"):
        network(reaction(glucose=-1, balanced_by=["heterotroph"]))


def test_closers_the_balances_cannot_tell_apart_are_refused():
    # Lactate and glucose both hold 4 electrons per carbon and no nitrogen.
    with pytest.raises(
        ConfigError, match="'glucose' holds carbon, nitrogen and electrons in the same proportions"
    ):
        network(reaction(ethanol=-1, balanced_by=["lactate", "glucose"]))


@pytest.mark.parametrize(
    ("process", "message"),
    [
        (reaction(glucose=-1, balanced_by=["glucose"]), "'glucose' already given a coefficient"),
        (
            reaction(glucose=-1, balanced_by=["oxygen", "carbon_dioxide", "ammonium", "nitrate"]),
            "only 3 balances",
        ),
        (reaction(glucos=-1), "'glucos' is not a component of this network; did you mean"),
        (reaction(glucose=-1, balanced_by=["oxygn"]), "did you mean 'oxygen'"),
        (reaction(), "is empty"),
        (reaction(glucose=-1, oxygen=0, balanced_by=["carbon_dioxide"]), "coefficient of zero"),
        (growth("glucose", "lactate", 1, []), "'glucose' is dissolved"),
        (growth("heterotroph", "heterotroph", 1, []), "both the biomass and its substrate"),
        (growth("heterotroph", "glucose", 0, []), "must be positive"),
        (growth("heterotroph", "glucose", -3.6, []), "must be positive"),
        (
            growth("heterotroph", "glucose", 3.6, [], products={"lactate": 0}),
            "lactate: must be positive",
        ),
        (
            growth("heterotroph", "glucose", 3.6, [], products={"glucose": 1}),
            "is this process's substrate",
        ),
    ],
)
def test_malformed_processes_are_refused(process, message):
    with pytest.raises(ConfigError, match=re.escape(message)):
        network(process)


# --- strict reading -----------------------------------------------------------


def with_first_process(**changes) -> dict:
    process = growth("heterotroph", "glucose", 3.6, CLOSE_AEROBICALLY, name="aerobic")
    process.update(changes)
    return raw_network(process)


@pytest.mark.parametrize(
    ("raw", "message"),
    [
        (
            with_first_process(yeild_mol_per_mol=3.6),
            "unknown field 'yeild_mol_per_mol'; did you mean 'yield_mol_per_mol'?",
        ),
        (
            with_first_process(yield_g_per_g=0.5),
            "unknown field 'yield_g_per_g'; the unit is part of the name, "
            "and this field is 'yield_mol_per_mol'",
        ),
        (with_first_process(rate_per_h=0.5), "unknown field 'rate_per_h'"),
        (with_first_process(yield_mol_per_mol=True), "expected a number, got true"),
        (with_first_process(yield_mol_per_mol="3.6"), "expected a number, got the text '3.6'"),
        (with_first_process(yield_mol_per_mol=float("nan")), "must be a finite number"),
        (with_first_process(kind="grow"), "expected 'growth' or 'reaction', got 'grow'"),
        (with_first_process(name="aerobic growth"), "is not a valid name"),
        (with_first_process(balanced_by="oxygen"), "expected a list of names"),
        (
            with_first_process(balanced_by=["oxygen", "oxygen"]),
            "oxygen listed more than once",
        ),
    ],
)
def test_mistakes_in_a_process_are_named(raw, message):
    with pytest.raises(ConfigError, match=re.escape(message)):
        network_from_dict(raw)


def test_a_process_without_its_kind_is_refused():
    raw = with_first_process()
    del raw["processes"][0]["kind"]
    with pytest.raises(ConfigError, match="expected 'growth' or 'reaction', got missing"):
        network_from_dict(raw)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"name": "../glucose"}, "is not a valid name"),
        ({"name": "1glucose"}, "is not a valid name"),
        ({"phase": "gas"}, "expected one of 'dissolved', 'particulate'"),
        ({"formula": "H2O"}, "Water and protons are never listed as components"),
        ({"formula": "NH4", "charge": 0}, "odd number of electrons"),
        ({"formule": "C6H12O6"}, "unknown field 'formule'; did you mean 'formula'?"),
    ],
)
def test_mistakes_in_a_component_are_named(change, message):
    raw = raw_network()
    raw["components"][0].update(change)
    with pytest.raises(ConfigError, match=re.escape(message)):
        network_from_dict(raw)


def test_names_must_be_unique():
    raw = raw_network()
    raw["components"].append(component("glucose"))
    with pytest.raises(ConfigError, match="names must be unique; glucose repeated"):
        network_from_dict(raw)
    twice = reaction("same", glucose=-1, balanced_by=["oxygen", "carbon_dioxide"])
    with pytest.raises(ConfigError, match="names must be unique; same repeated"):
        network(twice, twice)


def test_the_schema_version_is_required_and_checked():
    raw = raw_network()
    del raw["schema_version"]
    with pytest.raises(ConfigError, match="not a version 2 configuration"):
        network_from_dict(raw)
    with pytest.raises(ConfigError, match="3 is not supported"):
        network_from_dict(raw_network() | {"schema_version": 3})
    with pytest.raises(ConfigError, match="unknown field 'species'"):
        network_from_dict(raw_network() | {"species": []})


def test_a_version_one_ecosystem_file_is_pointed_at_marse_ecosystem():
    path = ROOT / "examples" / "experiments" / "two_species_ecosystem.json"
    with pytest.raises(ConfigError, match="marse ecosystem"):
        load_network(path)


def test_a_key_given_twice_in_a_file_is_refused(tmp_path):
    path = tmp_path / "twice.json"
    text = EXAMPLE.read_text("utf-8").replace(
        '"yield_mol_per_mol": 3.6,', '"yield_mol_per_mol": 3.6, "yield_mol_per_mol": 1.0,'
    )
    path.write_text(text, "utf-8")
    with pytest.raises(ConfigError, match="'yield_mol_per_mol' appears twice"):
        load_network(path)


def test_invalid_json_names_the_file_and_line(tmp_path):
    path = tmp_path / "broken.json"
    path.write_text('{\n  "schema_version": 2\n  "components": []\n}', "utf-8")
    with pytest.raises(ConfigError, match=r"broken\.json: invalid JSON .* at line 3"):
        load_network(path)


# --- properties ------------------------------------------------------------------


def random_network(seed: int) -> dict:
    """Random formulas and given rows, closed with oxygen, CO2 and ammonium.

    Those three always determine a row: their compositions (0, 0, -4),
    (1, 0, 0) and (0, 1, 0) span carbon, nitrogen and electrons.
    """
    rng = random.Random(seed)
    components = [component(name) for name in CLOSE_AEROBICALLY]
    for i in range(6):
        carbon, hydrogen = rng.randint(0, 6), rng.randint(0, 14)
        oxygen, nitrogen = rng.randint(0, 6), rng.randint(0, 2)
        charge = (hydrogen + nitrogen) % 2  # an even number of electrons
        if carbon == nitrogen == 0 and hydrogen - 2 * oxygen - charge == 0:
            carbon = 1  # something the balances can see: water and H+ are refused
        formula = "".join(
            f"{symbol}{count}"
            for symbol, count in (("C", carbon), ("H", hydrogen), ("N", nitrogen), ("O", oxygen))
            if count
        )
        phase = rng.choice(["dissolved", "particulate"])
        components.append({"name": f"s{i}", "phase": phase, "formula": formula, "charge": charge})
    processes = []
    for p in range(8):
        chosen = rng.sample([f"s{i}" for i in range(6)], rng.randint(1, 3))
        stoichiometry = {name: rng.choice([-1, 1]) * rng.randint(1, 40) / 10 for name in chosen}
        processes.append(reaction(f"p{p}", CLOSE_AEROBICALLY, **stoichiometry))
    return {"schema_version": 2, "components": components, "processes": processes}


@pytest.mark.parametrize("seed", range(40))
def test_random_networks_balance_exactly_down_to_the_atoms(seed):
    raw = random_network(seed)
    loaded = network_from_dict(raw)
    atoms = {}
    for item in raw["components"]:
        counts = dict.fromkeys("CHNO", 0)
        for symbol, digits in re.findall(r"([A-Z])(\d*)", item["formula"]):
            counts[symbol] += int(digits or 1)
        atoms[item["name"]] = (counts, item["charge"])
    for p in loaded.processes:
        # Recomputed from the atoms, independently of the loader's composition.
        for element in "CN":
            assert sum(c * atoms[n][0][element] for n, c in p.coefficients.items()) == 0
        oxygen = sum(c * atoms[n][0]["O"] for n, c in p.coefficients.items()) + p.water
        hydrogen = (
            sum(c * atoms[n][0]["H"] for n, c in p.coefficients.items()) + 2 * p.water + p.protons
        )
        charge = sum(c * atoms[n][1] for n, c in p.coefficients.items()) + p.protons
        assert (oxygen, hydrogen, charge) == (0, 0, 0)


@pytest.mark.parametrize("seed", range(10))
def test_the_float_matrices_balance_to_rounding(seed):
    # What the engine will integrate: each coefficient rounded once to float.
    loaded = network_from_dict(random_network(seed))
    stoichiometry, composition = loaded.stoichiometric_matrix(), loaded.composition_matrix()
    residual = np.abs(stoichiometry @ composition)
    scale = np.abs(stoichiometry) @ np.abs(composition)
    assert np.all(residual <= 16 * np.finfo(float).eps * scale)


@pytest.mark.invariance
def test_order_in_the_file_changes_nothing():
    raw = random_network(7)
    reference = {p.name: p.coefficients for p in network_from_dict(raw).processes}
    rng = random.Random(1)
    for _ in range(5):
        shuffled = json.loads(json.dumps(raw))
        rng.shuffle(shuffled["components"])
        rng.shuffle(shuffled["processes"])
        for process in shuffled["processes"]:
            rng.shuffle(process["balanced_by"])
        again = {p.name: p.coefficients for p in network_from_dict(shuffled).processes}
        assert again == reference


@pytest.mark.parametrize("seed", range(5))
def test_a_network_survives_the_trip_through_json(seed):
    loaded = network_from_dict(random_network(seed))
    assert network_from_dict(json.loads(json.dumps(loaded.to_dict()))) == loaded


def test_the_example_survives_the_trip_through_json_with_defaults_written_out():
    loaded = load_network(EXAMPLE)
    written = loaded.to_dict()
    assert written["components"][0]["charge"] == 0
    assert all("balanced_by" in p for p in written["processes"])
    assert network_from_dict(json.loads(json.dumps(written))) == loaded


# --- the schema documents itself ------------------------------------------------


def test_every_numeric_field_names_its_unit():
    for obj, fields in SCHEMA.items():
        for name, field in fields.items():
            if field.kind in NUMERIC_KINDS:
                assert name in DIMENSIONLESS or name.endswith(tuple(UNIT_SUFFIXES)), f"{obj}.{name}"


def test_the_documentation_lists_exactly_the_fields_the_reader_accepts():
    text = (ROOT / "docs" / "networks.md").read_text("utf-8")
    for obj, fields in SCHEMA.items():
        heading = f"### {obj.capitalize()} fields"
        assert heading in text, heading
        section = text.split(heading, 1)[1].split("\n#", 1)[0]
        documented = re.findall(r"^\| `([a-z_0-9]+)` \|", section, flags=re.MULTILINE)
        assert documented == list(fields), obj


# --- the example and the command line ---------------------------------------------


def test_the_example_derives_the_fermenters_lactate_from_glucose():
    loaded = load_network(EXAMPLE)
    fermentation = loaded.process("fermenter_growth_on_glucose")
    assert fermentation.coefficient("glucose") == F(-5, 3)
    assert fermentation.coefficient("lactate") == F(179, 60)
    assert loaded.process("heterotroph_growth_on_lactate").coefficient("lactate") == F(-2, 3)


def test_marse_check_prints_each_process_as_a_balanced_equation(capsys):
    assert main(["check", str(EXAMPLE)]) == 0
    out = capsys.readouterr().out
    assert "processes   4, each conserving carbon, nitrogen and electrons exactly" in out
    assert "0.277778 glucose + 0.616667 oxygen + 0.2 ammonium ->" in out
    assert "CH1.8O0.5N0.2" in out
    assert "NH4 (+1)" in out


def test_marse_check_refuses_an_unbalanced_network(tmp_path, capsys):
    raw = json.loads(EXAMPLE.read_text("utf-8"))
    respiration = raw["processes"][3]
    respiration["balanced_by"] = ["oxygen", "carbon_dioxide"]  # the nitrogen has nowhere to go
    path = tmp_path / "leaky.json"
    path.write_text(json.dumps(raw), "utf-8")
    assert main(["check", str(path)]) == 2
    out = capsys.readouterr().out
    assert "leaky.json is not a valid network" in out
    assert respiration["name"] in out
    assert "nothing in balanced_by carries nitrogen" in out


def test_marse_check_points_a_version_one_file_at_marse_ecosystem(capsys):
    experiment = ROOT / "examples" / "experiments" / "two_species_ecosystem.json"
    assert main(["check", str(experiment)]) == 2
    assert "marse ecosystem" in capsys.readouterr().out
