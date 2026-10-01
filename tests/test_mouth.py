"""The mouth over a column of plaque: clearance, the diet, conservation and the coupling.

The criteria set before Stage S1 was built that these check:

- G6a: with a constant flow and nothing taken up, the mouth's sugar falls by
  RESID / VMAX of the pool at every swallow, as Dawes's (1983) model says;
- G3 and G6b: the box and the mouth together conserve every quantity,
  counting what was secreted, eaten, swallowed and expelled, to rounding;
- G6c: the answer does not depend on how far the mouth runs ahead of the box,
  even while a rinse mixes the film fast.

The diet's rinses, drinks and foods each bring into the mouth what they state,
when they state it, and the food they leave on the teeth goes where it is put.
A food that is chewed adds the mouth's chewing flow while it lasts, and P7,
set before Stage S2 was built: sugar-free gum after a sugar rinse brings the
plaque's pH back sooner.
"""

import copy
import csv
import json
from pathlib import Path

import numpy as np
import pytest

import marse.core.reactive_transport as reactive_transport
from marse.chemistry import ChargeBalance
from marse.cli import main
from marse.core.config import ConfigError
from marse.oral import OralFluid, film_layers, renewal_per_h
from marse.schemas import experiment_from_dict, network_from_dict
from marse.spatial.grid import Grid

ROOT = Path(__file__).resolve().parents[1]
NETWORK = {
    "schema_version": 2,
    "pkw": 13.6,
    "components": [
        {"name": "sugar", "phase": "dissolved", "formula": "C6H12O6"},
        {
            "name": "lactate",
            "phase": "dissolved",
            "formula": "C3H6O3",
            "acid_base": {"pka": [3.86]},
        },
        {
            "name": "carbonate",
            "phase": "dissolved",
            "formula": "H2CO3",
            "acid_base": {"pka": [6.1, 10.0]},
        },
        {
            "name": "phosphate",
            "phase": "dissolved",
            "formula": "H3PO4",
            "acid_base": {"pka": [2.0, 6.8, 11.7]},
        },
        {"name": "potassium", "phase": "dissolved", "formula": "K", "charge": 1},
        {"name": "chloride", "phase": "dissolved", "formula": "Cl", "charge": -1},
        {
            "name": "carboxyl_groups",
            "phase": "particulate",
            "formula": "CH2O2",
            "acid_base": {"pka": [4.8]},
        },
        {"name": "bound_potassium", "phase": "particulate", "formula": "K", "charge": 1},
        {"name": "bacteria", "phase": "particulate", "formula": "CH1.8O0.5N0.2"},
    ],
    "processes": [
        {
            "name": "fermentation",
            "kind": "reaction",
            "stoichiometry_mol_per_mol": {"sugar": -1, "lactate": 2},
            "rate": {
                "maximum_per_h": 0.5,
                "proportional_to": "bacteria",
                "factors": [
                    {"component": "sugar", "form": "monod", "half_saturation_mol_per_m3": 1},
                    {"form": "ph", "ph_min": 4, "ph_optimum": 7, "ph_max": 9},
                ],
            },
        }
    ],
}
NAMES = [c["name"] for c in NETWORK["components"]]
BALANCE = ChargeBalance.of(network_from_dict(NETWORK))


def _neutral(composition, ph=7.0):
    c = np.array([composition.get(n, 0.0) for n in NAMES])
    net, _ = BALANCE.residual(c, np.array(1000.0 * 10**-ph))
    return {**composition, "potassium": composition.get("potassium", 0.0) - float(net)}


SALIVA = _neutral({"carbonate": 5.0, "phosphate": 5.0, "chloride": 20.0})
_KA, _H7 = 1000.0 * 10**-4.8, 1000.0 * 10**-7.0


def scene(sugar=0.0, *, duration_h=1 / 6, timestep_h=1 / 60, velocity=2.0, stimulated=2.0):
    """A column of 300 um of plaque under 100 um of film, renewed every 3 min at 2 mm/min."""
    raw = copy.deepcopy(NETWORK)
    colonies = [
        ("bacteria", 800.0),
        ("carboxyl_groups", 120.0),
        ("bound_potassium", 120.0 * _KA / (_KA + _H7)),
    ]
    raw.update(
        experiment_id="plaque_under_the_mouth",
        duration_h=duration_h,
        timestep_h=timestep_h,
        relative_tolerance=1e-3,
        absolute_tolerance_mol_per_m3=1e-6,
        initial_mol_per_m3=dict(SALIVA),
        domain={
            "voxels": [100],
            "voxel_um": 4,
            "diffusivity_m2_per_s": {
                "sugar": 2.3e-10,
                "lactate": 3.9e-10,
                "carbonate": 1.1e-9,
                "phosphate": 4.0e-10,
                "potassium": 1.1e-9,
                "chloride": 1.1e-9,
            },
            "colonies": [
                {
                    "component": name,
                    "center_um": [],
                    "radius_um": 300,
                    "concentration_mol_per_m3": amount,
                }
                for name, amount in colonies
            ],
            "film": {
                "thickness_um": 100,
                "velocity_mm_per_min": velocity,
                "plaque_length_mm": 6,
            },
            "mouth": {
                "saliva_mol_per_m3": dict(SALIVA),
                "resting_volume_ml": 0.77,
                "swallow_volume_ml": 1.07,
                "unstimulated_flow_ml_per_min": 0.3,
                "stimulated_flow_ml_per_min": stimulated,
                "stimulus": "sugar",
                "stimulus_half_mol_per_m3": 50,
                "plaque_area_cm2": 2,
                "initial_mol_per_m3": {**SALIVA, "sugar": sugar},
            },
        },
    )
    return raw


# --- the film ------------------------------------------------------------------------


def test_the_film_is_renewed_fastest_at_its_surface_and_not_at_all_below_it():
    config = experiment_from_dict(scene())
    rate = renewal_per_h(config.domain.grid, config.domain.film)
    assert np.all(rate[:75] == 0.0)
    assert np.all(np.diff(rate[75:]) > 0)
    mean = rate[75:].mean()  # u_bar / l on average: 2 mm/min over 6 mm is 20 per hour
    assert mean == pytest.approx(20.0, rel=2e-3)


# --- G6a: clearance ------------------------------------------------------------------


def test_g6a_the_mouth_clears_by_the_same_share_at_every_swallow():
    """Constant flow, a film that hardly exchanges: sugar falls by H_resid / H_vmax per swallow."""
    raw = scene(sugar=100.0, duration_h=10 / 60, timestep_h=1 / 60, velocity=1e-9, stimulated=0.0)
    del raw["domain"]["mouth"]["stimulus"]
    del raw["domain"]["mouth"]["stimulus_half_mol_per_m3"]
    raw["processes"] = []
    config = experiment_from_dict(raw)
    result = reactive_transport.run(config)
    fluid = OralFluid(config.domain.mouth, config.domain.film, config.network.component_names)
    share = fluid.thickness_um(fluid.resting_m3) / fluid.thickness_um(fluid.full_m3)
    sugar = result.mouth["sugar_mol_per_m3"]
    assert result.mouth["swallows"][-1] == 10  # one a minute at 0.3 mL/min
    expected = 100.0 * share ** np.arange(11)
    np.testing.assert_allclose(sugar, expected, rtol=1e-9)


# --- G3 and G6b: conservation --------------------------------------------------------


def test_g3_the_box_and_the_mouth_conserve_what_was_secreted_and_swallowed():
    result = reactive_transport.run(experiment_from_dict(scene(sugar=584.0)))
    mouth = result.manifest.outputs["mouth"]
    assert mouth["swallows"] > 5
    for quantity, entry in mouth["balance"].items():
        assert entry["largest_relative_residual"] < 1e-12, quantity
    for quantity, entry in result.manifest.outputs["balance"].items():
        assert entry["largest_relative_residual"] < 1e-12, quantity
    assert result.ph["substratum_min"][-1] < 7.0  # sugar reached the plaque and made acid
    assert mouth["final_mol_per_m3"]["lactate"] > 0.0  # and acid reached the mouth


def test_a_column_at_rest_under_resting_saliva_stays_as_it_is():
    result = reactive_transport.run(experiment_from_dict(scene()))
    np.testing.assert_allclose(result.ph["box_min"], 7.0, atol=1e-9)
    np.testing.assert_allclose(result.ph["box_max"], 7.0, atol=1e-9)
    np.testing.assert_allclose(result.mouth["ph"], 7.0, atol=1e-9)
    assert result.mouth["swallows"][-1] == 10


# --- G6c: how far the mouth runs ahead ----------------------------------------------------


def test_g6c_the_answer_does_not_depend_on_how_far_the_mouth_runs_ahead(monkeypatch):
    raw = scene(sugar=584.0, duration_h=6 / 60, timestep_h=0.5 / 60)
    long = reactive_transport.run(experiment_from_dict(raw))
    monkeypatch.setattr(reactive_transport, "_LONGEST_SPAN_S", 5.0)
    short = reactive_transport.run(experiment_from_dict(raw))
    assert short.manifest.outputs["mouth"]["spans"] > 2 * long.manifest.outputs["mouth"]["spans"]
    assert np.max(np.abs(long.ph["substratum_mean"] - short.ph["substratum_mean"])) < 1e-3
    np.testing.assert_allclose(
        long.mouth["sugar_mol_per_m3"], short.mouth["sugar_mol_per_m3"], rtol=1e-3
    )


# --- the command line, the schema, replay ------------------------------------------------


def test_a_run_under_the_mouth_writes_the_mouth_and_replays(tmp_path, capsys):
    path = tmp_path / "scene.json"
    path.write_text(json.dumps(scene(sugar=584.0, duration_h=3 / 60)), encoding="utf-8")
    assert main(["run", str(path), "-o", str(tmp_path / "out")]) == 0
    with (tmp_path / "out" / "mouth.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert list(rows[0])[:4] == ["time_h", "volume_ml", "flow_ml_per_min", "swallows"]
    assert float(rows[0]["volume_ml"]) == pytest.approx(0.77)
    assert float(rows[0]["sugar_mol_per_m3"]) == pytest.approx(584.0)
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["models"]["mouth"] == reactive_transport.MOUTH_VERSION
    capsys.readouterr()
    assert (
        main(["replay", str(tmp_path / "out" / "manifest.json"), "-o", str(tmp_path / "again")])
        == 0
    )
    assert "reproduced the recorded results exactly" in capsys.readouterr().out


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda d: d.pop("mouth"), "film and the mouth that renews it come together; mouth"),
        (lambda d: d["film"].update(thickness_um=102), "not a whole number of voxels"),
        (lambda d: d["film"].update(thickness_um=400), "the film fills the box"),
        (lambda d: d["mouth"].update(swallow_volume_ml=0.5), "more than resting_volume_ml"),
        (lambda d: d["mouth"].pop("stimulus"), "needs the stimulus"),
        (lambda d: d["mouth"].update(stimulus="bacteria"), "is particulate"),
        (lambda d: d["mouth"].update(plaque_area_cm2=100), "no less than the resting volume"),
        (lambda d: d.update(bulk_mol_per_m3={"sugar": 1}), "the liquid comes from the mouth"),
        (lambda d: d["mouth"]["saliva_mol_per_m3"].update(bacteria=1), "is particulate"),
    ],
)
def test_impossible_mouths_are_refused(change, message):
    raw = scene()
    change(raw["domain"])
    with pytest.raises(ConfigError, match=message):
        experiment_from_dict(raw)


def test_a_mouth_is_written_back_as_read():
    config = experiment_from_dict(scene(sugar=10.0))
    again = experiment_from_dict(json.loads(json.dumps(config.to_dict())))
    assert again == config
    assert config.to_dict()["domain"]["mouth"]["initial_mol_per_m3"]["sugar"] == 10.0


def test_the_grid_of_a_scene_is_what_the_film_assumes():
    assert Grid((100,), 4.0).size_um[-1] == 400.0


# --- the diet ----------------------------------------------------------------------------

AREA_M2 = 2e-4  # the scene's plaque, 2 cm2
FOOD = {"name": "food_sugar", "phase": "particulate", "formula": "C6H12O6"}
DISSOLVING = {
    "name": "dissolving",
    "kind": "reaction",
    "stoichiometry_mol_per_mol": {"food_sugar": -1, "sugar": 1},
    "rate": {"maximum_per_h": 6, "proportional_to": "food_sugar"},
}
RINSE = {
    "kind": "rinse",
    "start_h": 0.0,
    "duration_min": 1,
    "volume_ml": 10,
    "composition_mol_per_m3": {"sugar": 584},
}


def dieted(diet, *, pocket=False, **kwargs):
    """The scene with a diet; with ``pocket``, food that sticks to the teeth and dissolves."""
    raw = scene(**kwargs)
    if pocket:
        raw["components"].append(dict(FOOD))
        raw["processes"].append(copy.deepcopy(DISSOLVING))
    raw["domain"]["diet"] = copy.deepcopy(diet)
    return raw


def test_g3_every_intake_is_booked_as_eaten_and_conserved():
    diet = [
        RINSE,
        {
            "kind": "drink",
            "start_h": 2 / 60,
            "duration_min": 2,
            "volume_ml": 20,
            "composition_mol_per_m3": {"sugar": 300, "potassium": 10, "chloride": 10},
        },
        {
            "kind": "food",
            "start_h": 4.5 / 60,
            "duration_min": 2.5,
            "released_mmol": {"sugar": 5},
            "retained": {"component": "food_sugar", "amount_mol_per_m2": 0.05},
        },
    ]
    result = reactive_transport.run(
        experiment_from_dict(dieted(diet, pocket=True, duration_h=9 / 60, timestep_h=0.5 / 60))
    )
    mouth = result.manifest.outputs["mouth"]
    assert mouth["intakes"] == 3
    assert result.manifest.models["diet"] == reactive_transport.DIET_VERSION
    eaten = mouth["eaten_mol_per_m2"]
    sugar = (584 * 10e-6 + 300 * 20e-6 + 5e-3) / AREA_M2
    assert eaten["sugar"] == pytest.approx(sugar, rel=1e-12)
    assert eaten["potassium"] == pytest.approx(10 * 20e-6 / AREA_M2, rel=1e-12)
    assert eaten["food_sugar"] == pytest.approx(0.05, rel=1e-12)
    assert mouth["expelled_mol_per_m2"]["sugar"] > 0.5 * 584 * 10e-6 / AREA_M2
    assert mouth["swallowed_mol_per_m2"]["sugar"] > 0.0
    for quantity, entry in mouth["balance"].items():
        assert entry["largest_relative_residual"] < 1e-12, quantity
    for quantity, entry in result.manifest.outputs["balance"].items():
        assert entry["largest_relative_residual"] < 1e-12, quantity
    # The box's imports count the pocket, placed in it, as well as what crossed into the film.
    assert result.imported_mol_per_m2[-1][result.component_names.index("food_sugar")] == (
        pytest.approx(0.05, rel=1e-12)
    )


def test_a_rinse_is_held_without_swallowing_then_expelled_to_the_resting_volume():
    rinse = dict(RINSE, start_h=0.3 / 60)  # from 0.3 to 1.3 minutes, between records
    result = reactive_transport.run(
        experiment_from_dict(dieted([rinse], duration_h=3 / 60, timestep_h=0.25 / 60))
    )
    t = np.round(result.times_h * 60, 9)
    volume, swallows = result.mouth["volume_ml"], result.mouth["swallows"]
    held = (t >= 0.5) & (t <= 1.25)
    assert np.all(volume[held] > 10.77)
    assert np.all(np.diff(volume[held]) > 0)  # the glands keep secreting into it
    assert np.all(swallows[held] == swallows[t == 0.25])
    assert np.all(volume[t >= 1.5] <= 1.07)
    assert result.manifest.outputs["mouth"]["intakes"] == 1


def test_expelling_a_rinse_keeps_the_resting_share_and_changes_no_concentration():
    config = experiment_from_dict(scene())
    fluid = OralFluid(config.domain.mouth, config.domain.film, config.network.component_names)
    fluid.take(10e-6)
    held = fluid.thickness_um(fluid.volume_m3)
    kept = fluid.expel()
    assert fluid.volume_m3 == fluid.resting_m3
    assert kept == pytest.approx(fluid.thickness_um(fluid.resting_m3) / held, rel=1e-15)


def test_a_drink_is_swallowed_as_it_fills_the_mouth_and_sets_its_sugar():
    drink = {
        "kind": "drink",
        "start_h": 0.0,
        "duration_min": 2,
        "volume_ml": 20,
        "composition_mol_per_m3": {"sugar": 300},
    }
    result = reactive_transport.run(
        experiment_from_dict(dieted([drink], duration_h=2 / 60, timestep_h=0.5 / 60))
    )
    # 10 mL a minute and the saliva with it: a swallow for every 0.3 mL.
    assert result.mouth["swallows"][-1] > 2 * 10 / 0.3
    # The mouth holds the drink diluted by the stimulated flow: c = q c_drink / (q + Q).
    flow = result.mouth["flow_ml_per_min"][-1]
    expected = 10 * 300 / (10 + flow)
    assert result.mouth["sugar_mol_per_m3"][-1] == pytest.approx(expected, rel=0.02)


def test_a_sweet_releases_its_sugar_steadily_without_liquid():
    sweet = {"kind": "food", "start_h": 0.0, "duration_min": 5, "released_mmol": {"sugar": 10}}
    result = reactive_transport.run(
        experiment_from_dict(dieted([sweet], duration_h=5 / 60, timestep_h=0.5 / 60))
    )
    mouth = result.manifest.outputs["mouth"]
    assert mouth["eaten_mol_per_m2"]["sugar"] == pytest.approx(10e-3 / AREA_M2, rel=1e-12)
    assert mouth["eaten_mol_per_m2"]["potassium"] == 0.0
    # Tasted, the sugar raises the flow well above the resting 0.3 mL per minute.
    assert np.all(result.mouth["flow_ml_per_min"][2:] > 1.0)
    assert np.all(np.diff(result.mouth["sugar_mol_per_m3"][:3]) > 0)


def test_food_left_on_the_teeth_goes_into_the_film_over_its_region():
    raw = dieted(
        [
            {
                "kind": "food",
                "start_h": 0.0,
                "duration_min": 0.5,
                "released_mmol": {"sugar": 1},
                "retained": {
                    "component": "food_sugar",
                    "amount_mol_per_m2": 0.2,
                    "region_um": [0, 8],
                },
            }
        ],
        pocket=True,
        duration_h=1 / 60,
        timestep_h=0.5 / 60,
    )
    raw["domain"]["voxels"] = [4, 100]  # four columns of 4 um: the region covers two
    for colony in raw["domain"]["colonies"]:
        colony["center_um"] = [8]  # 300 um of plaque under all four
    config = experiment_from_dict(raw)
    result = reactive_transport.run(config)
    food = result.final_state[result.component_names.index("food_sugar")]
    layers = film_layers(config.domain.grid, config.domain.film)
    assert np.all(food[:2, -layers:] > 0.0)
    assert np.all(food[2:] == 0.0)
    assert np.all(food[:, :-layers] == 0.0)
    eaten = result.manifest.outputs["mouth"]["eaten_mol_per_m2"]["food_sugar"]
    assert eaten == pytest.approx(0.2 * 8 / 16, rel=1e-12)  # per m2 of the whole substratum
    # Half an hour's first-order release at 6 per hour has dissolved some of it.
    sugar = result.manifest.outputs["totals_mol_per_m2"]
    assert sugar["food_sugar"] < eaten


def test_mixing_during_an_intake_brings_its_sugar_to_the_plaque_faster():
    def plaque_sugar(mixing):
        rinse = dict(RINSE, mixing_per_s=mixing)
        result = reactive_transport.run(
            experiment_from_dict(dieted([rinse], duration_h=1 / 60, timestep_h=0.5 / 60))
        )
        return result.final_state[NAMES.index("sugar"), :75].mean()

    assert plaque_sugar(1.0) > 1.5 * plaque_sugar(0.0)


def test_a_diet_that_starts_after_the_run_changes_nothing():
    plain = reactive_transport.run(experiment_from_dict(scene(sugar=100.0, duration_h=2 / 60)))
    later = dieted([dict(RINSE, start_h=1.0)], sugar=100.0, duration_h=2 / 60)
    dieted_run = reactive_transport.run(experiment_from_dict(later))
    digest = plain.manifest.outputs["final_state_sha256"]
    assert dieted_run.manifest.outputs["final_state_sha256"] == digest
    assert dieted_run.manifest.outputs["mouth"]["intakes"] == 0


def test_g6c_a_rinse_that_mixes_the_film_fast_does_not_depend_on_the_spans(monkeypatch):
    raw = dieted([RINSE], duration_h=6 / 60, timestep_h=0.5 / 60)
    long = reactive_transport.run(experiment_from_dict(raw))
    monkeypatch.setattr(reactive_transport, "_LONGEST_SPAN_S", 5.0)
    short = reactive_transport.run(experiment_from_dict(raw))
    assert short.manifest.outputs["mouth"]["spans"] > 2 * long.manifest.outputs["mouth"]["spans"]
    assert np.max(np.abs(long.ph["substratum_mean"] - short.ph["substratum_mean"])) < 1e-3
    np.testing.assert_allclose(
        long.mouth["sugar_mol_per_m3"], short.mouth["sugar_mol_per_m3"], rtol=1e-3
    )


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda d: d.pop("mouth"), "a diet is taken into the mouth"),
        (lambda d: d["diet"][0].update(kind="snack"), "kind"),
        (lambda d: d["diet"][0].pop("volume_ml"), "a rinse needs volume_ml"),
        (lambda d: d["diet"][0].update(released_mmol={"sugar": 1}), "give it as composition"),
        (lambda d: d["diet"][0].update(start_h=-1), "must not be negative"),
        (lambda d: d["diet"][0].update(duration_min=0), "must be positive"),
        (lambda d: d["diet"][0].update(mixing_per_s=-1), "must not be negative"),
        (lambda d: d["diet"][0]["composition_mol_per_m3"].update(bacteria=1), "a rinse holds"),
        (lambda d: d["diet"][1].update(volume_ml=5), "food adds no liquid"),
        (lambda d: d["diet"][1].pop("released_mmol"), "food needs released_mmol"),
        (lambda d: d["diet"][1].update(start_h=0.5 / 60), "before intake 0 ends"),
        (
            lambda d: d["diet"][1]["retained"].update(component="sugar"),
            "food left on the teeth is particulate",
        ),
        (lambda d: d["diet"][1]["retained"].update(component="bacteria"), "no process consumes"),
        (lambda d: d["diet"][1]["retained"].update(region_um=[0, 4]), "needs 0 bounds"),
        (lambda d: d["diet"][1]["retained"].update(amount_mol_per_m2=0), "must be positive"),
        (lambda d: d["diet"][0].update(chewing=True), "only food is chewed"),
        (lambda d: d["diet"][1].update(chewing=True), "give the mouth chewing_flow_ml_per_min"),
        (lambda d: d["diet"][1].update(chewing="yes"), "expected true or false"),
        (lambda d: d["mouth"].update(chewing_flow_ml_per_min=-1), "must not be negative"),
    ],
)
def test_impossible_diets_are_refused(change, message):
    diet = [
        RINSE,
        {
            "kind": "food",
            "start_h": 2 / 60,
            "duration_min": 1,
            "released_mmol": {"sugar": 1},
            "retained": {"component": "food_sugar", "amount_mol_per_m2": 0.01},
        },
    ]
    raw = dieted(diet, pocket=True)
    change(raw["domain"])
    with pytest.raises(ConfigError, match=message):
        experiment_from_dict(raw)


def test_a_diet_is_written_back_as_read_and_replays(tmp_path, capsys):
    diet = [
        RINSE,
        {
            "kind": "food",
            "start_h": 1.5 / 60,
            "duration_min": 1,
            "released_mmol": {"sugar": 2},
            "mixing_per_s": 0.5,
            "retained": {"component": "food_sugar", "amount_mol_per_m2": 0.01},
        },
    ]
    raw = dieted(diet, pocket=True, duration_h=3 / 60)
    config = experiment_from_dict(raw)
    again = experiment_from_dict(json.loads(json.dumps(config.to_dict())))
    assert again == config
    written = config.to_dict()["domain"]["diet"]
    assert written[0]["mixing_per_s"] == 1.0  # the default, written out
    assert written[1]["retained"]["region_um"] == []
    path = tmp_path / "diet.json"
    path.write_text(json.dumps(raw), encoding="utf-8")
    assert main(["check", str(path)]) == 0
    out = " ".join(capsys.readouterr().out.split())
    assert "rinse at 0 min for 1 min: 10 mL holding sugar 584 mol per m3" in out
    assert "food at 1.5 min for 1 min: releases sugar 2 mmol; mixes the film at 0.5 per s" in out
    assert "leaves 0.01 mol per m2 of food_sugar on the teeth" in out
    assert main(["run", str(path), "-o", str(tmp_path / "out")]) == 0
    out = " ".join(capsys.readouterr().out.split())
    assert "counting what was secreted, eaten, swallowed and expelled" in out
    assert "counting what crossed the film (food left on the teeth included)" in out
    assert (
        main(["replay", str(tmp_path / "out" / "manifest.json"), "-o", str(tmp_path / "again")])
        == 0
    )
    assert "reproduced the recorded results exactly" in capsys.readouterr().out


# --- chewing ----------------------------------------------------------------------------

GUM = {"kind": "food", "start_h": 0.0, "duration_min": 10, "released_mmol": {}, "chewing": True}
STEPHAN_RINSE = ROOT / "examples" / "environments" / "oral" / "stephan_rinse.json"


def chewed(diet, *, chewing_flow=1.0, **kwargs):
    """The scene with a diet, a mouth that chews, and stimulated saliva richer in carbonate."""
    raw = dieted(diet, **kwargs)
    mouth = raw["domain"]["mouth"]
    mouth["chewing_flow_ml_per_min"] = chewing_flow
    mouth["stimulated_saliva_mol_per_m3"] = _neutral(
        {"carbonate": 15.0, "phosphate": 4.0, "chloride": 20.0}
    )
    return raw


def test_chewing_adds_its_flow_while_it_lasts_and_the_saliva_comes_stimulated():
    raw = chewed([GUM], duration_h=15 / 60, timestep_h=1 / 60)
    result = reactive_transport.run(experiment_from_dict(raw))
    minutes = result.times_h * 60
    flow = result.mouth["flow_ml_per_min"]
    np.testing.assert_allclose(flow[(minutes > 0.5) & (minutes < 9.5)], 1.3, rtol=1e-12)
    np.testing.assert_allclose(flow[minutes > 10.5], 0.3, rtol=1e-12)
    # At 1.3 mL a minute the mouth swallows about four times a minute, against once at rest.
    swallows = result.mouth["swallows"]
    assert swallows[10] >= 40
    assert swallows[-1] - swallows[10] <= 6
    # While chewing, the glands secrete saliva half way to stimulated: more carbonate.
    carbonate = result.mouth["carbonate_mol_per_m3"]
    assert carbonate[9] > carbonate[0] + 2.0
    outputs = result.manifest.outputs
    for balance in (outputs["balance"], outputs["mouth"]["balance"]):
        for quantity, entry in balance.items():
            assert entry["largest_relative_residual"] < 1e-12, quantity


def test_chewing_is_written_back_as_read_and_described(tmp_path, capsys):
    raw = chewed([GUM], duration_h=2 / 60)
    config = experiment_from_dict(raw)
    written = config.to_dict()["domain"]
    assert written["diet"][0]["chewing"] is True
    assert written["mouth"]["chewing_flow_ml_per_min"] == 1.0
    assert experiment_from_dict(json.loads(json.dumps(config.to_dict()))) == config
    # A mouth that does not chew writes neither field, so S1's configurations read as before.
    plain = experiment_from_dict(scene()).to_dict()["domain"]
    assert "chewing_flow_ml_per_min" not in plain["mouth"]
    path = tmp_path / "gum.json"
    path.write_text(json.dumps(raw), encoding="utf-8")
    assert main(["check", str(path)]) == 0
    out = " ".join(capsys.readouterr().out.split())
    assert "1 more while chewing" in out
    assert "food at 0 min for 10 min: releases nothing mmol, chewed" in out


@pytest.mark.slow
def test_p7_gum_chewed_after_a_sugar_rinse_brings_the_plaque_back_sooner():
    raw = json.loads(STEPHAN_RINSE.read_text("utf-8"))
    raw["duration_h"] = 0.25
    raw["domain"]["mouth"]["chewing_flow_ml_per_min"] = 1.0
    plain = reactive_transport.run(experiment_from_dict(raw))
    gum = dict(GUM, start_h=2 / 60, duration_min=20)
    raw["domain"]["diet"].append(gum)
    chewing = reactive_transport.run(experiment_from_dict(raw))
    rinse, gummed = plain.ph["substratum_mean"], chewing.ph["substratum_mean"]
    # The rinse alone is still falling at 15 minutes; gum has stopped the fall and turned it.
    assert rinse[-1] < 5.0
    assert gummed[-1] > 6.0
    assert gummed.min() > rinse.min() + 0.5
