"""Plaque that spreads in a column, wears at its surface, and is brushed off.

The criteria set before Stage S2 was built, checked here (docs/validation.md,
"Plaque that spreads"):

- P2: the solid never overfills a voxel, and its front is sharp: full voxels,
  then at most one partly filled;
- P3: a film that grows at a constant specific rate mu and wears at u follows
  dL/dt = mu L - u (Wanner and Gujer 1986);
- P4: a brushing removes exactly its share, from the surface down.

And conservation: what detaches, wears away or is brushed off is booked, so
both ledgers close to rounding.
"""

import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest

import marse.core.reactive_transport as reactive_transport
from marse.biofilm.spreading import SPREADING_VERSION, Spreading, remap
from marse.cli import main
from marse.core.config import ConfigError
from marse.schemas import experiment_from_dict

ROOT = Path(__file__).resolve().parents[1]
RINSE = ROOT / "examples" / "environments" / "oral" / "stephan_rinse.json"
POCKET = ROOT / "examples" / "environments" / "oral" / "pocket.json"
PACKING = 800.0  # mol/m3


def _reference(column, occupying, packing, moving, keep):
    """The remap done slowly, voxel by voxel: what each old voxel's solid becomes."""
    phi = (column[occupying] / packing[:, None]).sum(axis=0)
    top = np.cumsum(phi)
    bottom = top - phi
    cut = min(keep, top[-1], column.shape[1])
    new = column.copy()
    new[moving] = 0.0
    detached = np.zeros(column.shape[0])
    for k in range(column.shape[1]):
        if phi[k] == 0:
            new[moving, k] += column[moving, k]
            continue
        for j in range(column.shape[1]):
            low, high = max(bottom[k], j), min(top[k], j + 1, cut)
            if high > low:
                new[moving, j] += column[moving, k] * (high - low) / phi[k]
        above = max(0.0, top[k] - max(bottom[k], cut))
        detached[moving] += column[moving, k] * above / phi[k]
    return new, detached


def test_the_remap_packs_the_solid_conserves_it_and_cuts_it_exactly():
    rng = np.random.default_rng(2)
    occupying, packing = np.array([0, 1]), np.array([800.0, 400.0])
    moving = np.array([0, 1, 2])  # two that fill space, one carried
    for _ in range(500):
        voxels = int(rng.integers(3, 40))
        column = np.zeros((4, voxels))
        filled = int(rng.integers(1, voxels))
        column[0, :filled] = rng.uniform(0, 900, filled) * (rng.random(filled) < 0.9)
        column[1, :filled] = rng.uniform(0, 300, filled) * (rng.random(filled) < 0.9)
        column[2, :filled] = rng.uniform(0, 160, filled)  # left where no solid carries it
        column[3] = rng.uniform(0, 5, voxels)  # dissolved: never moves
        height = float((column[occupying] / packing[:, None]).sum())
        keep = float(rng.choice([voxels, height * rng.uniform(0, 1.2), rng.uniform(0, voxels)]))
        new, detached = remap(column, occupying, packing, moving, keep)
        expected, gone = _reference(column, occupying, packing, moving, keep)
        scale = column[:3].sum()
        np.testing.assert_allclose(new, expected, rtol=0, atol=1e-12 * scale)
        np.testing.assert_allclose(detached, gone, rtol=0, atol=1e-12 * scale)
        assert np.all(new >= 0)
        assert np.all(detached >= 0)
        assert np.array_equal(new[3], column[3])
        residual = new[:3].sum(axis=1) + detached[:3] - column[:3].sum(axis=1)
        assert np.all(np.abs(residual) <= 1e-15 * scale)
        phi = (new[occupying] / packing[:, None]).sum(axis=0)
        assert np.all(phi <= 1 + 1e-12)  # P2
        partial = np.flatnonzero((phi > 1e-12) & (phi < 1 - 1e-12))
        assert partial.size <= 1
        assert float(phi.sum()) == pytest.approx(min(keep, height, voxels), abs=1e-9)


# --- a film in the laboratory: growth, wear and the maximum height ---------------------


def _film(L0_um, *, mu=0.1, wear=4.0, maximum=200.0, hours=5.0, records_h=0.5, voxel=2.5, box=None):
    voxels = round((box or maximum) / voxel)
    return {
        "schema_version": 2,
        "experiment_id": "spreading-film",
        "components": [
            {"name": "nutrient", "phase": "dissolved", "formula": "CH1.8O0.5N0.2"},
            {"name": "biomass", "phase": "particulate", "formula": "CH1.8O0.5N0.2"},
            {"name": "walls", "phase": "particulate", "formula": "CH2O2"},
        ],
        "processes": [
            {
                "name": "growth",
                "kind": "reaction",
                "stoichiometry_mol_per_mol": {"nutrient": -1, "biomass": 1},
                "rate": {
                    "maximum_per_h": mu,
                    "proportional_to": "biomass",
                    "assumed_in_excess": ["nutrient"],
                },
            }
        ],
        "initial_mol_per_m3": {"nutrient": 1e6},
        "duration_h": hours,
        "timestep_h": records_h,
        "relative_tolerance": 1e-6,
        "absolute_tolerance_mol_per_m3": 1e-12,
        "domain": {
            "voxels": [voxels],
            "voxel_um": voxel,
            "bulk_mol_per_m3": {"nutrient": 1e6},
            "diffusivity_m2_per_s": {"nutrient": 1e-9},
            "colonies": [
                {
                    "component": "biomass",
                    "center_um": [],
                    "radius_um": L0_um,
                    "concentration_mol_per_m3": PACKING,
                },
                {
                    "component": "walls",
                    "center_um": [],
                    "radius_um": L0_um,
                    "concentration_mol_per_m3": 160.0,
                },
            ],
            "plaque": {
                "packing_mol_per_m3": {"biomass": PACKING},
                "carried": ["walls"],
                "maximum_um": maximum,
                "wear_um_per_h": wear,
            },
        },
    }


@pytest.mark.parametrize(("L0", "hours"), [(20.0, 5.0), (60.0, 10.0)])
def test_p3_a_film_growing_against_wear_follows_wanner_and_gujer(L0, hours):
    mu, wear = 0.1, 4.0
    result = reactive_transport.run(experiment_from_dict(_film(L0, mu=mu, wear=wear, hours=hours)))
    L_star = wear / mu  # unstable: below it the film wears away, above it the film grows
    exact = L_star + (L0 - L_star) * np.exp(mu * result.times_h)
    np.testing.assert_allclose(result.plaque["thickness_um"], exact, rtol=2e-6)
    outputs = result.manifest.outputs
    assert result.manifest.models["spreading"] == SPREADING_VERSION
    assert outputs["plaque"]["spreading"] == SPREADING_VERSION
    for quantity, entry in outputs["balance"].items():
        assert entry["largest_relative_residual"] < 1e-13, quantity


def test_wear_alone_takes_the_surface_off_at_its_velocity():
    raw = _film(30.0, mu=0.0, wear=5.0, hours=8.0, records_h=1.0)
    result = reactive_transport.run(experiment_from_dict(raw))
    expected = np.maximum(30.0 - 5.0 * result.times_h, 0.0)
    np.testing.assert_allclose(result.plaque["thickness_um"], expected, atol=1e-9)
    worn = result.manifest.outputs["plaque"]["detached_mol_per_m2"]
    assert worn["biomass"] == pytest.approx(PACKING * 30e-6, rel=1e-12)  # all of it, in mol/m2
    assert worn["walls"] == pytest.approx(160.0 * 30e-6, rel=1e-12)  # carried away with it


def test_growth_past_the_maximum_height_is_detached_and_booked():
    raw = _film(40.0, mu=0.2, wear=0.0, maximum=60.0, hours=6.0, records_h=1.0)
    result = reactive_transport.run(experiment_from_dict(raw))
    thickness = result.plaque["thickness_um"]
    assert thickness.max() == pytest.approx(60.0, abs=1e-9)
    assert thickness[-1] == pytest.approx(60.0, abs=1e-9)
    # The column reaches 60 um at t = ln(1.5) / mu; from then on, what the full 60 um grows
    # is detached. Cut after every step, the excess grows a little first: first order.
    reached = math.log(60.0 / 40.0) / 0.2
    detached = result.manifest.outputs["plaque"]["detached_mol_per_m2"]["biomass"]
    continuous = 0.2 * PACKING * 60e-6 * (6.0 - reached)
    assert detached > continuous
    assert detached == pytest.approx(continuous, rel=2e-3)
    for quantity, entry in result.manifest.outputs["balance"].items():
        assert entry["largest_relative_residual"] < 1e-13, quantity


def test_carried_components_move_with_the_solid():
    raw = _film(40.0, mu=0.2, wear=0.0, maximum=200.0, hours=2.0, records_h=1.0)
    result = reactive_transport.run(experiment_from_dict(raw))
    names = result.component_names
    biomass, walls = (result.final_state[names.index(n)] for n in ("biomass", "walls"))
    held = biomass > 0
    # Growth makes biomass but no walls, so every voxel thins its walls by the same factor:
    # the share it holds depends only on the solid's age where it sits, never on a jump.
    assert np.all(walls[~held] == 0)
    assert walls.sum() == pytest.approx(160.0 * 40 / 2.5, rel=1e-12)  # all still there


# --- brushing under the mouth --------------------------------------------------------------


def _brushed(fraction=None, kind="brushing", at_h=0.25):
    raw = json.loads(RINSE.read_text("utf-8"))
    raw["duration_h"] = 0.5
    raw["timestep_h"] = 0.05
    raw["domain"]["diet"] = []
    raw["domain"]["plaque"] = {
        "packing_mol_per_m3": {"bacteria": 800.0},
        "carried": ["carboxyl_groups", "bound_potassium"],
    }
    event = {"kind": kind, "start_h": at_h}
    if fraction is not None:
        event["removes_fraction"] = fraction
    raw["domain"]["hygiene"] = [event]
    return raw


def test_p4_a_brushing_takes_its_share_off_from_the_surface_and_the_mouth_expels_it():
    result = reactive_transport.run(experiment_from_dict(_brushed()))
    thickness = result.plaque["thickness_um"]
    times = result.times_h
    assert thickness[times < 0.25 - 1e-9] == pytest.approx(150.0, abs=1e-9)
    assert thickness[times > 0.25 + 1e-9] == pytest.approx(150.0 * (1 - 0.42), abs=1e-9)
    outputs = result.manifest.outputs
    removed = outputs["plaque"]["removed_mol_per_m2"]
    assert removed["bacteria"] == pytest.approx(0.42 * 800.0 * 150e-6, rel=1e-12)
    assert removed["carboxyl_groups"] == pytest.approx(0.42 * 160.0 * 150e-6, rel=1e-12)
    mouth = outputs["mouth"]
    assert mouth["cleanings"] == 1
    assert mouth["expelled_mol_per_m2"]["bacteria"] == removed["bacteria"]
    for balance in (outputs["balance"], mouth["balance"]):
        for quantity, entry in balance.items():
            assert entry["largest_relative_residual"] < 1e-13, quantity


def test_a_brushing_takes_the_same_share_of_food_left_on_the_teeth():
    raw = json.loads(POCKET.read_text("utf-8"))
    raw["duration_h"] = 0.2
    raw["domain"]["plaque"] = {
        "packing_mol_per_m3": {"bacteria": 800.0},
        "carried": ["carboxyl_groups", "bound_potassium"],
    }
    raw["domain"]["hygiene"] = [{"kind": "brushing", "start_h": 0.1}]
    result = reactive_transport.run(experiment_from_dict(raw))
    food = result.component_names.index("food_sugar")
    at = int(np.flatnonzero(np.isclose(result.times_h, 0.1))[0])  # recorded before the brush
    before = result.totals_mol_per_m2[at][food]
    assert before > 0.01  # most of the 0.02 mol/m2 is still on the teeth at 6 minutes
    outputs = result.manifest.outputs
    assert outputs["mouth"]["removed_mol_per_m2"]["food_sugar"] == pytest.approx(
        0.42 * before, rel=1e-12
    )
    for balance in (outputs["balance"], outputs["mouth"]["balance"]):
        for quantity, entry in balance.items():
            assert entry["largest_relative_residual"] < 1e-13, quantity


def test_a_flossing_removes_the_share_it_states():
    result = reactive_transport.run(experiment_from_dict(_brushed(0.25, kind="flossing")))
    assert result.plaque["thickness_um"][-1] == pytest.approx(150.0 * 0.75, abs=1e-9)


def test_a_plaque_that_does_not_spread_runs_as_before():
    raw = json.loads(RINSE.read_text("utf-8"))
    raw["duration_h"] = 0.05
    plain = reactive_transport.run(experiment_from_dict(raw))
    raw["domain"]["plaque"] = {
        "packing_mol_per_m3": {"bacteria": 800.0},
        "carried": ["carboxyl_groups", "bound_potassium"],
    }
    packed = reactive_transport.run(experiment_from_dict(raw))
    # S1's rinse grows nothing, so its plaque stays packed where it was: the same run.
    np.testing.assert_allclose(packed.final_state, plain.final_state, rtol=1e-12, atol=1e-15)
    assert plain.plaque is None
    assert "plaque" not in plain.manifest.outputs
    assert "spreading" not in plain.manifest.models


# --- the schema ------------------------------------------------------------------------


def _two_dimensional(raw):
    raw["domain"]["voxels"] = [4, 80]
    for colony in raw["domain"]["colonies"]:
        colony["center_um"] = [5.0]


def _refused(change, match):
    raw = _film(30.0)
    change(raw)
    with pytest.raises(ConfigError, match=match):
        experiment_from_dict(raw)


@pytest.mark.parametrize(
    ("change", "match"),
    [
        (lambda r: _two_dimensional(r), "in a column only"),
        (lambda r: r["domain"]["plaque"].update(packing_mol_per_m3={"nutrient": 1.0}), "dissolved"),
        (lambda r: r["domain"]["plaque"].update(packing_mol_per_m3={}), "at least one"),
        (lambda r: r["domain"]["plaque"].update(packing_mol_per_m3={"biomass": 0}), "positive"),
        (lambda r: r["domain"]["plaque"].update(carried=["biomass"]), "already fills"),
        (lambda r: r["domain"]["plaque"].update(carried=["walls", "walls"]), "more than once"),
        (lambda r: r["domain"]["plaque"].update(carried=["sugar"]), "sugar"),
        (lambda r: r["domain"]["plaque"].update(maximum_um=500), "fit below the top of the box"),
        (lambda r: r["domain"]["plaque"].update(wear_um_per_h=-1), "must not be negative"),
        (lambda r: r["domain"]["plaque"].update(wear_um_per_s=1), "wear_um_per_h"),
        (lambda r: r["domain"].update(hygiene=[{"kind": "brushing", "start_h": 1}]), "a mouth"),
    ],
)
def test_impossible_plaques_are_refused(change, match):
    _refused(change, match)


@pytest.mark.parametrize(
    ("event", "match"),
    [
        ({"kind": "flossing", "start_h": 0.1}, "only brushing has a measured default"),
        ({"kind": "brushing", "start_h": 0.1, "removes_fraction": 0}, "more than 0"),
        ({"kind": "brushing", "start_h": 0.1, "removes_fraction": 1.5}, "at most 1"),
        ({"kind": "brushing", "start_h": -1}, "must not be negative"),
        ({"kind": "scraping", "start_h": 0.1}, "scraping"),
    ],
)
def test_impossible_cleanings_are_refused(event, match):
    raw = _brushed()
    raw["domain"]["hygiene"] = [event]
    with pytest.raises(ConfigError, match=match):
        experiment_from_dict(raw)


def test_cleanings_come_in_time_order_and_need_a_plaque():
    raw = _brushed()
    raw["domain"]["hygiene"] = [
        {"kind": "brushing", "start_h": 0.3},
        {"kind": "brushing", "start_h": 0.1},
    ]
    with pytest.raises(ConfigError, match="time order"):
        experiment_from_dict(raw)
    raw = _brushed()
    raw["domain"].pop("plaque")
    with pytest.raises(ConfigError, match="give the domain a plaque"):
        experiment_from_dict(raw)
    raw = _brushed()
    raw["domain"]["plaque"]["maximum_um"] = 200  # into the film
    with pytest.raises(ConfigError, match="fit below the film, 150 um up"):
        experiment_from_dict(raw)


def test_food_left_on_the_teeth_is_not_part_of_the_plaque():
    raw = json.loads(POCKET.read_text("utf-8"))
    raw["domain"]["plaque"] = {"packing_mol_per_m3": {"bacteria": 800.0}, "carried": ["food_sugar"]}
    with pytest.raises(ConfigError, match="'food_sugar' is food left on the teeth"):
        experiment_from_dict(raw)


def test_an_initial_plaque_taller_than_its_maximum_is_refused():
    raw = _film(30.0, maximum=20.0, box=100.0)
    with pytest.raises(ConfigError, match="30 um of solid, more than the plaque's maximum of 20"):
        reactive_transport.run(experiment_from_dict(raw))


def test_a_plaque_and_its_cleanings_read_back_as_written():
    raw = _brushed()
    config = experiment_from_dict(raw)
    written = config.domain.to_dict()
    assert written["plaque"]["maximum_um"] == pytest.approx(150.0)
    assert written["hygiene"] == [{"kind": "brushing", "start_h": 0.25, "removes_fraction": 0.42}]
    again = copy.deepcopy(raw)
    again["domain"] = written
    assert experiment_from_dict(again).domain == config.domain
    plain = json.loads(RINSE.read_text("utf-8"))
    assert "plaque" not in experiment_from_dict(plain).domain.to_dict()


def test_a_run_whose_plaque_spreads_replays_bit_for_bit(tmp_path, capsys):
    path = tmp_path / "film.json"
    path.write_text(json.dumps(_film(20.0, hours=1.0)), "utf-8")
    out = tmp_path / "run"
    assert main(["run", str(path), "-o", str(out)]) == 0
    assert (out / "plaque.csv").read_text("utf-8").startswith("time_h,thickness_um,")
    assert main(["replay", str(out / "manifest.json")]) == 0
    assert "reproduced the recorded results exactly" in capsys.readouterr().out


def test_check_and_run_describe_the_plaque_and_its_cleanings(tmp_path, capsys):
    raw = _brushed()
    raw["domain"]["plaque"]["wear_um_per_h"] = 2.0
    path = tmp_path / "brushed.json"
    path.write_text(json.dumps(raw), "utf-8")
    assert main(["check", str(path)]) == 0
    out = capsys.readouterr().out
    assert "plaque      spreads up the column (displacement_1d_v1): filled by bacteria" in out
    assert "at most 150 um high, detached above; wears 2 um per hour at its surface" in out
    assert "brushing    at 15 min: takes 42% of the plaque off from its surface down" in out
    assert main(["run", str(path), "-o", str(tmp_path / "run")]) == 0
    out = capsys.readouterr().out
    assert "plaque      86.21 um at the end" in out  # (150 - 0.5) x 0.58 - 0.5
    assert "taken off by 1 cleaning(s), mol per m2: bacteria 0.05023" in out


def test_spreading_needs_a_column():
    spreading = Spreading(
        occupying=np.array([0]),
        packing_mol_per_m3=np.array([PACKING]),
        moving=np.array([0]),
        maximum_voxels=10.0,
        voxel_um=2.5,
    )
    box = np.zeros((1, 10))
    box[0, :4] = PACKING
    assert spreading.height_um(box) == pytest.approx(10.0)
    removed, taken = spreading.remove(box, 0.5)
    assert spreading.height_um(removed) == pytest.approx(5.0)
    assert taken[0] == pytest.approx(2 * PACKING)
