"""Surfaces, part two: cells in the liquid binding to the substratum, in the engine.

The engine is checked against the closed-form kinetics in one, two and three
dimensions, against its own Langmuir and jamming limits, patch by patch
against single columns, against the ledger, and against finite differences
of its Jacobian. The schema refuses every impossible scene, and the two
scenes run from the command line. docs/validation.md, "Adhesion to surfaces".
"""

import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest

from marse.cli import main
from marse.core.config import ConfigError
from marse.core.ledger import Ledger
from marse.core.reactive_transport import build_model, surface_exchange
from marse.core.reactive_transport import run as run_in_space
from marse.microbes.adhesion import JAMMING_COVERAGE
from marse.schemas.experiment import experiment_from_dict
from marse.schemas.formula import QUANTITIES
from marse.validation.analytical import adhesion_kinetics

ROOT = Path(__file__).resolve().parents[1]
LAB = ROOT / "examples" / "environments" / "lab" / "flow_chamber.json"
DENTAL = ROOT / "examples" / "environments" / "dental" / "dental_surfaces.json"
CELL = {"name": "", "phase": "particulate", "formula": "CH1.8O0.5N0.2"}


def scene(voxels=(4,), patches=None, adhesion=None, duration_h=4.0, **overrides) -> dict:
    """One species binding in buffer: nothing dissolved, nothing grows."""
    lateral = len(voxels) - 1
    size = [5.0 * n for n in voxels[:lateral]]
    whole = [x for s in size for x in (0.0, s)]
    raw = {
        "schema_version": 2,
        "experiment_id": "binding",
        "components": [dict(CELL, name="bound"), dict(CELL, name="cells")],
        "processes": [],
        "duration_h": duration_h,
        "timestep_h": 0.25,
        "relative_tolerance": 1e-6,
        "absolute_tolerance_mol_per_m3": 1e-3,
        "domain": {
            "voxels": list(voxels),
            "voxel_um": 5.0,
            "diffusivity_m2_per_s": {},
            "substratum": {
                "conditioning_film": "none",
                "patches": patches or [{"material": "glass", "region_um": whole}],
            },
            "liquid": {"temperature_c": 37.0, "viscosity_mpa_s": 0.692},
            "flow": {"wall_shear_rate_per_s": 15.0, "distance_from_inlet_mm": 20.0},
            "suspension": [
                {
                    "reversible": "bound",
                    "attached": "cells",
                    "cells_per_ml": 3e8,
                    "cell_diameter_um": 0.9,
                    "carbon_fmol_per_cell": 7.0,
                    "blocked_area_um2": 10.0,
                }
            ],
            "adhesion": adhesion
            or [
                {
                    "attached": "cells",
                    "material": "glass",
                    "efficiency": 0.5,
                    "detachment_per_h": 30.0,
                    "locking_per_h": 90.0,
                }
            ],
        },
    }
    raw.update(overrides)
    return raw


def _rates(config, material="glass"):
    """j0 and n_J, per cm2, from the configuration, as the closed form takes them."""
    exchange = surface_exchange(config)
    species = exchange.species[0]
    face = config.domain.substratum.mask(material)
    efficiency = float(np.asarray(species.efficiency)[face].ravel()[0])
    arrival = efficiency * species.arrival_um_per_h * species.cells_per_um3 * 1e8  # per cm2 per h
    jammed = JAMMING_COVERAGE / species.blocked_area_um2 * 1e8
    return arrival, jammed


# --- the kinetics, in the engine -------------------------------------------------------


@pytest.mark.parametrize("voxels", [(4,), (3, 4), (2, 2, 4)])
def test_binding_follows_the_closed_form_kinetics_in_every_dimension(voxels):
    config = experiment_from_dict(scene(voxels))
    result = run_in_space(config)
    arrival, jammed = _rates(config)
    reversible, locked = adhesion_kinetics(result.times_h, arrival, jammed, 30.0, 90.0)
    got = result.surface["glass_cells_cells_per_cm2"]
    np.testing.assert_allclose(got, reversible + locked, rtol=1e-5, atol=1e-6 * jammed)


def test_without_locking_binding_settles_on_the_langmuir_balance():
    rule = {"attached": "cells", "material": "glass", "efficiency": 0.5}
    raw = scene(adhesion=[dict(rule, detachment_per_h=30.0, locking_per_h=0.0)], duration_h=6.0)
    config = experiment_from_dict(raw)
    result = run_in_space(config)
    arrival, jammed = _rates(config)
    settled = arrival / (arrival / jammed + 30.0)  # n = j0 / (j0 / n_J + k_off)
    assert result.surface["glass_cells_cells_per_cm2"][-1] == pytest.approx(settled, rel=1e-6)


def test_coverage_reaches_the_jamming_limit_and_never_passes_it():
    rule = {"attached": "cells", "material": "glass", "efficiency": 1.0}
    raw = scene(adhesion=[dict(rule, detachment_per_h=0.0, locking_per_h=200.0)], duration_h=48.0)
    raw["timestep_h"] = 2.0
    config = experiment_from_dict(raw)
    exchange = surface_exchange(config)
    model = build_model(config)
    names = config.network.component_names
    state = config.domain.initial_state(names, config.initial_mol_per_m3, 0)
    blocked = []
    for _ in range(24):
        state, *_ = model.integrate(state, 2.0, relative_tolerance=1e-6, absolute_tolerance=1e-3)
        blocked.append(float(exchange.coverage(state).max()))
    assert max(blocked) <= JAMMING_COVERAGE * (1 + 1e-9)
    assert blocked[-1] == pytest.approx(JAMMING_COVERAGE, rel=1e-3)


def test_each_patch_evolves_as_its_own_column():
    two = [
        {"material": "glass", "region_um": [0.0, 10.0]},
        {"material": "coated", "region_um": [10.0, 20.0]},
    ]
    rules = [
        {
            "attached": "cells",
            "material": "glass",
            "efficiency": 0.5,
            "detachment_per_h": 30.0,
            "locking_per_h": 90.0,
        },
        {
            "attached": "cells",
            "material": "coated",
            "efficiency": 0.3,
            "detachment_per_h": 10.0,
            "locking_per_h": 60.0,
        },
    ]
    slice_ = build_model(experiment_from_dict(scene((4, 4), patches=two, adhesion=rules)))
    y = np.zeros((2, 4, 4))
    columns = {}
    for rule in rules:
        column_rule = dict(rule, material="glass")
        columns[rule["material"]] = (
            build_model(experiment_from_dict(scene((4,), adhesion=[column_rule]))),
            np.zeros((2, 4)),
        )
    for _ in range(12):  # the same fixed steps everywhere, so only the geometry could differ
        y, *_ = slice_.step(y, 0.1, np.abs(y), 1e-3, 1e-6)
        for name, (model, state) in columns.items():
            columns[name] = (model, model.step(state, 0.1, np.abs(state), 1e-3, 1e-6)[0])
    for i, name in ((0, "glass"), (1, "glass"), (2, "coated"), (3, "coated")):
        np.testing.assert_allclose(y[:, i, :], columns[name][1], rtol=1e-12, atol=1e-12)


# --- conservation and positivity -------------------------------------------------------


def test_the_ledger_counts_what_binds_and_what_detaches():
    config = experiment_from_dict(scene((3, 4)))
    result = run_in_space(config)
    for quantity, entry in result.manifest.outputs["balance"].items():
        assert entry["largest_relative_residual"] < 1e-12, quantity
        assert entry["imported"] > 0, quantity  # the cells came in through the substratum
    # Nothing grows here, so every cell on the surface crossed the substratum.
    totals = result.totals_mol_per_m2[-1]
    imported = result.imported_mol_per_m2[-1]
    assert totals.sum() == pytest.approx(imported.sum(), rel=1e-12)


def test_hour_long_steps_stay_positive_and_nothing_is_clipped():
    config = experiment_from_dict(scene((3, 4)))
    model = build_model(config)
    names = config.network.component_names
    state = config.domain.initial_state(names, config.initial_mol_per_m3, 0)
    ledger = Ledger(config.network.composition_matrix(), state, QUANTITIES, tolerance=1e-12)
    for step in range(1, 7):
        state, imports, *_ = model.step(state, 1.0, np.abs(state), 1e-3, 1e-4)
        assert np.all(state >= 0)
        ledger.exchange(imports)
        ledger.check(state, step=step, time_h=float(step))


def test_the_exchange_jacobian_matches_finite_differences():
    config = experiment_from_dict(json.loads(DENTAL.read_text("utf-8")))
    exchange = surface_exchange(config)
    rng = np.random.default_rng(1)
    names = config.network.component_names
    c = rng.uniform(0.0, 50.0, size=(len(names), 16, 4, 8))
    c[rng.random(c.shape) < 0.2] *= -1  # undershoots, which nothing responds to
    height = config.domain.grid.voxel_um
    jacobian = exchange.exchange_jacobian(c)
    assert exchange.coverage(np.abs(c)).max() > JAMMING_COVERAGE  # faces on both sides of jamming
    for k in range(len(names)):
        d = np.zeros_like(c)
        d[k, ..., 0] = 1e-6 * np.maximum(1.0, np.abs(c[k, ..., 0]))
        numerical = (exchange.exchange(c + d) - exchange.exchange(c - d)) / (2 * d[k, ..., 0])
        np.testing.assert_allclose(
            jacobian[:, k], numerical / height, rtol=1e-6, atol=1e-9 * np.abs(jacobian).max()
        )
    locking = exchange.locking_jacobian(c)
    for k in range(len(names)):
        d = np.zeros_like(c)
        d[k] = 1e-6 * np.maximum(1.0, np.abs(c[k]))
        numerical = (exchange.locking(c + d) - exchange.locking(c - d)) / (2 * d[k])
        np.testing.assert_allclose(locking[:, k], numerical, rtol=1e-6, atol=1e-6)


# --- the schema --------------------------------------------------------------------------


def test_a_domain_without_a_surface_writes_what_it_always_wrote():
    raw = scene()
    for key in ("substratum", "liquid", "flow", "suspension", "adhesion"):
        del raw["domain"][key]
    written = experiment_from_dict(raw).domain.to_dict()
    assert set(written) == {
        "voxels",
        "voxel_um",
        "bulk_mol_per_m3",
        "diffusivity_m2_per_s",
        "colonies",
        "random_colonies",
    }  # so a run without a surface keeps its run id and its checksums
    with_surface = experiment_from_dict(scene())
    assert experiment_from_dict(json.loads(json.dumps(with_surface.to_dict()))) == with_surface


def _edited(change):
    raw = scene()
    change(raw)
    return raw


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda r: r["domain"].pop("flow"), "flow is missing"),
        (
            lambda r: r["components"][0].update(formula="C6H12O6"),
            "must have the same formula",
        ),
        (
            lambda r: (
                r["components"][0].update(phase="dissolved"),
                r["domain"]["diffusivity_m2_per_s"].update(bound=1e-9),
            ),
            "is dissolved",
        ),
        (lambda r: r["domain"]["suspension"][0].update(attached="bound"), "must be two"),
        (
            lambda r: r["domain"]["adhesion"][0].update(material="steel"),
            "not a material of the substratum",
        ),
        (lambda r: r["domain"]["adhesion"][0].update(efficiency=1.5), "between 0 and 1"),
        (lambda r: r["domain"]["adhesion"][0].update(detachment_per_h=-1), "must not be negative"),
        (
            lambda r: r["domain"]["adhesion"].append(copy.deepcopy(r["domain"]["adhesion"][0])),
            "stated twice",
        ),
        (lambda r: r["domain"].update(adhesion=[]), "missing 'cells' on 'glass'"),
        (lambda r: r["domain"]["liquid"].update(temperature_c=-300), "above absolute zero"),
        (lambda r: r["domain"]["flow"].update(wall_shear_rate_per_s=0), "must be positive"),
        (
            lambda r: r["domain"]["substratum"]["patches"].append(
                {"material": "steel", "region_um": []}
            ),
            "overlaps",
        ),
        (lambda r: r["domain"]["suspension"][0].update(cell_diameter_um=0), "must be positive"),
    ],
)
def test_an_impossible_scene_is_refused_with_the_reason(change, message):
    with pytest.raises(ConfigError, match=message):
        experiment_from_dict(_edited(change))


def test_reversibly_bound_cells_take_part_in_no_process():
    raw = json.loads(DENTAL.read_text("utf-8"))
    raw["processes"][0]["biomass"] = "s_oralis_reversible"
    # Growing makes it a species, which needs a role before the scene is read.
    reversible = next(c for c in raw["components"] if c["name"] == "s_oralis_reversible")
    reversible["oxygen_role"] = "aerotolerant"
    with pytest.raises(ConfigError, match="reversibly bound cells only detach or lock"):
        experiment_from_dict(raw)


# --- the scenes --------------------------------------------------------------------------


def test_the_lab_scene_binds_at_a_rate_inside_the_published_range():
    config = experiment_from_dict(json.loads(LAB.read_text("utf-8")))
    for material in ("glass", "saliva_coated_glass"):
        arrival, _ = _rates(config, material)
        assert 0 < arrival / 3600 < 2.9e3  # per cm2 per s (Sjollema, Busscher, Weerkamp 1988)


def test_the_dental_scene_binds_less_to_zirconia_by_the_calibrated_ratio():
    raw = json.loads(DENTAL.read_text("utf-8"))
    raw["duration_h"] = 2.0
    raw["domain"]["voxels"] = [8, 2, 4]  # the same 80 x 20 x 40 um, in coarser voxels
    raw["domain"]["voxel_um"] = 10.0
    raw["absolute_tolerance_mol_per_m3"] = 1e-3  # the ratio does not hang on resolving traces
    result = run_in_space(experiment_from_dict(raw))
    for species in ("s_oralis", "s_sanguinis"):
        titanium = result.surface[f"titanium_{species}_cells_per_cm2"][-1]
        zirconia = result.surface[f"zirconia_{species}_cells_per_cm2"][-1]
        assert zirconia / titanium == pytest.approx(0.63, rel=1e-2)
        assert result.surface[f"enamel_{species}_cells_per_cm2"][-1] == pytest.approx(titanium)


@pytest.mark.slow
def test_a_day_on_the_dental_surfaces_covers_titanium_more_than_zirconia():
    result = run_in_space(experiment_from_dict(json.loads(DENTAL.read_text("utf-8"))))
    covered = {m: result.surface[f"{m}_covered"][-1] for m in ("titanium", "zirconia")}
    # 19.3% against 12.1% in the mouth (Scarano et al. 2004): a ratio of 0.63
    assert 0.6 < covered["zirconia"] / covered["titanium"] < 0.72
    for entry in result.manifest.outputs["balance"].values():
        assert entry["largest_relative_residual"] < 1e-12


# --- the command line ----------------------------------------------------------------------


def test_marse_run_writes_the_surface_record_and_replays(tmp_path, capsys):
    output = tmp_path / "lab"
    assert main(["run", str(LAB), "-o", str(output)]) == 0
    out = capsys.readouterr().out
    assert "bound cells per cm2" in out
    assert "saliva_coated_glass" in out
    header = (output / "surface.csv").read_text("utf-8").splitlines()[0].split(",")
    assert header == [
        "time_h",
        "glass_s_oralis_cells_per_cm2",
        "saliva_coated_glass_s_oralis_cells_per_cm2",
        "glass_covered",
        "saliva_coated_glass_covered",
    ]
    assert main(["replay", str(output / "manifest.json")]) == 0
    assert "reproduced the recorded results exactly" in capsys.readouterr().out


def test_marse_check_previews_binding_on_each_material(capsys):
    assert main(["check", str(DENTAL)]) == 0
    out = capsys.readouterr().out
    assert "substratum: enamel 25%, titanium 25%, zirconia 25%, pmma 25%" in out
    assert "on zirconia: binds" in out
    assert math.isfinite(float(out.split("delivered at ")[1].split()[0]))
