"""The mouth over a column of plaque: clearance, conservation and the coupling.

The criteria set before Stage S1 was built that these check:

- G6a: with a constant flow and nothing taken up, the mouth's sugar falls by
  RESID / VMAX of the pool at every swallow, as Dawes's (1983) model says;
- G3 and G6b: the box and the mouth together conserve every quantity,
  counting what was secreted and swallowed, to rounding;
- G6c: the answer does not depend on how far the mouth runs ahead of the box.
"""

import copy
import csv
import json

import numpy as np
import pytest

import marse.core.reactive_transport as reactive_transport
from marse.chemistry import ChargeBalance
from marse.cli import main
from marse.core.config import ConfigError
from marse.oral import OralFluid, renewal_per_h
from marse.schemas import experiment_from_dict, network_from_dict
from marse.spatial.grid import Grid

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
