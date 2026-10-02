"""Reactions and transport together, in space: the spatial engine of Stage 2c.

The scheme is checked against its own exact algebra on a diffusion mode, the
engine against the well-mixed engine and against itself across dimensions,
and every run against its ledger, which now counts what crosses the top face.
docs/validation.md, "Reactions and transport in space".
"""

import copy
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from marse.cli import main
from marse.core.config import ConfigError
from marse.core.implicit import GAMMA, ReactionTransport
from marse.core.ledger import Ledger
from marse.core.provenance import Manifest
from marse.core.reactive_transport import build_model
from marse.core.reactive_transport import run as run_in_space
from marse.core.simulation import ConservationError
from marse.core.well_mixed import run as run_well_mixed
from marse.microbes.kinetics import compile_rates, process_rates, rate_jacobian
from marse.schemas.experiment import (
    ReactiveTransportConfig,
    WellMixedConfig,
    experiment_from_dict,
)
from marse.schemas.formula import QUANTITIES
from marse.spatial.grid import Grid
from marse.spatial.transport import Diffusion
from marse.validation.analytical import cosine_mode_rate

ROOT = Path(__file__).resolve().parents[1]
NETWORK = ROOT / "examples" / "networks" / "glucose_cross_feeding.json"
BULK = {"glucose": 5.0, "oxygen": 0.21, "ammonium": 5.0, "carbon_dioxide": 1.0}
# Effective diffusivities in a biofilm at 37 C (docs/parameters.md), m2/s.
DIFFUSIVITY = {
    "glucose": 2.3e-10,
    "oxygen": 1.13e-9,
    "ammonium": 1.47e-9,
    "carbon_dioxide": 1.10e-9,
    "lactate": 3.9e-10,
}


def spatial(voxels=(24,), voxel_um=4.0, duration_h=0.25, colonies=None, **overrides) -> dict:
    raw = json.loads(NETWORK.read_text("utf-8"))
    raw["experiment_id"] = "test-in-space"
    raw["relative_tolerance"] = 1e-3  # accuracy is tested where it matters, speed elsewhere
    raw["initial_mol_per_m3"] = dict(BULK)
    raw["duration_h"] = duration_h
    raw["timestep_h"] = min(0.05, duration_h)
    raw.pop("record_interval_h", None)
    lateral = len(voxels) - 1
    raw["domain"] = {
        "voxels": list(voxels),
        "voxel_um": voxel_um,
        "bulk_mol_per_m3": dict(BULK),
        "diffusivity_m2_per_s": dict(DIFFUSIVITY),
        "colonies": colonies
        if colonies is not None
        else [
            {
                "component": "heterotroph",
                "center_um": [voxel_um * voxels[a] / 2 for a in range(lateral)],
                "radius_um": 40.0,
                "concentration_mol_per_m3": 800.0,
            },
            {
                "component": "fermenter",
                "center_um": [voxel_um * voxels[a] / 2 for a in range(lateral)],
                "radius_um": 24.0,
                "concentration_mol_per_m3": 400.0,
            },
        ],
    }
    raw.update(overrides)
    return raw


SCATTERED = [
    {"component": "fermenter", "count": 2, "radius_um": 8.0, "concentration_mol_per_m3": 1.0}
]


def refused(raw: dict, message: str) -> None:
    with pytest.raises(ConfigError, match=message):
        experiment_from_dict(raw)


# --- the domain in the version 2 schema -----------------------------------------------


def test_a_domain_makes_a_document_run_in_space_and_round_trips():
    config = experiment_from_dict(spatial())
    assert isinstance(config, ReactiveTransportConfig)
    assert config.kind == "reactive_transport"
    assert config.relative_tolerance == 1e-3
    default = spatial()
    del default["relative_tolerance"]
    assert experiment_from_dict(default).relative_tolerance == 1e-4  # the default in space
    again = experiment_from_dict(json.loads(json.dumps(config.to_dict())))
    assert again == config
    assert config.to_dict()["domain"]["bulk_mol_per_m3"]["lactate"] == 0.0  # zeros written out
    without = spatial()
    del without["domain"]
    assert isinstance(experiment_from_dict(without), WellMixedConfig)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda d: d["diffusivity_m2_per_s"].pop("lactate"), "'lactate' has none"),
        (
            lambda d: d["diffusivity_m2_per_s"].update(heterotroph=1e-10),
            "'heterotroph' is particulate and does not diffuse",
        ),
        (lambda d: d["bulk_mol_per_m3"].update(fermenter=1.0), "'fermenter' is particulate"),
        (lambda d: d["colonies"][0].update(component="glucose"), "'glucose' is dissolved"),
        (lambda d: d["colonies"][0].update(radius_um=0.5), "covers no voxel"),
        (lambda d: d.update(voxels=[8, 8, 8, 8]), "one, two or three axes"),
        (lambda d: d.update(voxels=[8, 1]), "at least two voxels"),
        (lambda d: d.update(voxel_um=0), "must be positive"),
        (lambda d: d.update(voxel_mm=4), "this field is 'voxel_um', in micrometres"),
        (lambda d: d.update(random_colonies=SCATTERED), "a one-dimensional column has no surface"),
    ],
)
def test_domains_that_describe_no_possible_run_are_refused(change, message):
    raw = spatial()
    change(raw["domain"])
    refused(raw, message)


def test_a_colony_outside_the_box_is_refused():
    raw = spatial(voxels=(8, 16))
    raw["domain"]["colonies"][0]["center_um"] = [40.0]
    refused(raw, "lies outside the box")


def test_random_colonies_are_placed_by_the_seed_and_replaced_by_it():
    group = [
        {"component": "fermenter", "count": 3, "radius_um": 6.0, "concentration_mol_per_m3": 5.0}
    ]
    raw = spatial(voxels=(8, 8, 8), colonies=[])
    raw["domain"]["random_colonies"] = group
    first = experiment_from_dict(raw)
    names = first.network.component_names
    state = first.domain.initial_state(names, first.initial_mol_per_m3, first.seed)
    again = first.domain.initial_state(names, first.initial_mol_per_m3, first.seed)
    np.testing.assert_array_equal(state, again)
    other = first.domain.initial_state(names, first.initial_mol_per_m3, seed=7)
    assert not np.array_equal(state, other)
    assert first.domain.random_streams == ("colonies",)
    assert np.count_nonzero(state[names.index("fermenter")]) > 0


# --- the scheme's own algebra ---------------------------------------------------------


def diffusion_only(grid: Grid, diffusivity: float, bulk: float) -> ReactionTransport:
    return ReactionTransport(
        Diffusion(grid, np.array([diffusivity]), np.array([bulk])),
        np.zeros((1, 1)),
        lambda c: np.zeros((1, *c.shape[1:])),
        lambda c: np.zeros((1, 1, *c.shape[1:])),
    )


@pytest.mark.numerical
def test_a_diffusion_mode_decays_exactly_as_the_scheme_says_and_converges_at_second_order():
    # 2 + mode, with the bulk held at 2: the offset is an exact steady state, so
    # the mode evolves alone, and the scheme must multiply it by exactly
    # R(z) = (1 + (1 - 2 gamma) z) / (1 - gamma z)^2 per step, z = h * rate.
    grid = Grid((8, 6, 12), 3.0)
    diffusivity = 2.0e5
    k = [2 * np.pi / 24.0, 2 * np.pi / 18.0, np.pi / (2 * 36.0)]
    x, y, z = (grid.centers_um(a) for a in range(3))
    mode = (
        np.cos(k[0] * x)[:, None, None]
        * np.cos(k[1] * y)[None, :, None]
        * np.cos(k[2] * z)[None, None, :]
    )
    rate = cosine_mode_rate(np.array(k), diffusivity, grid.voxel_um)
    model = diffusion_only(grid, diffusivity, 2.0)
    span = 2.0 / abs(rate)  # two e-foldings: the non-stiff regime, where the order shows
    errors = []
    for n in (4, 8, 16, 32):
        h = span / n
        state = 2.0 + mode[None]
        for _ in range(n):
            state, *_ = model.step(state, h, np.abs(state), 1e-14, 1e-12)
        amplitude = np.sum((state[0] - 2.0) * mode) / np.sum(mode * mode)
        z_ = h * rate
        stability = (1 + (1 - 2 * GAMMA) * z_) / (1 - GAMMA * z_) ** 2
        assert amplitude == pytest.approx(stability**n, rel=1e-9)
        errors.append(abs(amplitude - math.exp(rate * span)))
    orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert np.all(orders > 1.9), orders


# --- the engine -----------------------------------------------------------------------


def test_without_diffusion_every_voxel_is_a_well_mixed_box():
    raw = spatial(voxels=(2, 4), duration_h=0.5)
    raw["domain"]["diffusivity_m2_per_s"] = dict.fromkeys(DIFFUSIVITY, 0.0)
    raw["relative_tolerance"] = 1e-6
    raw["domain"]["colonies"] = [
        {
            "component": "heterotroph",
            "center_um": [4.0],
            "radius_um": 20.0,
            "concentration_mol_per_m3": 50.0,
        }
    ]
    config = experiment_from_dict(raw)
    result = run_in_space(config)
    names = config.network.component_names
    state = config.domain.initial_state(names, config.initial_mol_per_m3, 0)
    for voxel in np.ndindex(*config.domain.grid.shape):
        start = {n: float(state[(j, *voxel)]) for j, n in enumerate(names)}
        box = copy.deepcopy(raw)
        del box["domain"]
        box["initial_mol_per_m3"] = start
        box["relative_tolerance"] = 1e-6
        mixed = run_well_mixed(experiment_from_dict(box))
        for j, n in enumerate(names):
            assert result.final_state[(j, *voxel)] == pytest.approx(
                mixed.final_mol_per_m3[n], rel=1e-4, abs=1e-7
            ), (voxel, n)


def test_a_laterally_uniform_box_behaves_as_one_column():
    column = experiment_from_dict(spatial(voxels=(16,)))
    box = experiment_from_dict(spatial(voxels=(2, 2, 16), colonies=[]))
    names = column.network.component_names
    column_model, box_model = build_model(column), build_model(box)
    one = column.domain.initial_state(names, column.initial_mol_per_m3, 0)  # layers, in 1-D
    one, *_ = column_model.integrate(one, 0.02, relative_tolerance=1e-4, absolute_tolerance=1e-9)
    box_state = np.broadcast_to(one[:, None, None, :], (len(names), 2, 2, 16)).copy()
    for _ in range(4):  # the same steps in both, so any difference is the discretisation's
        one, *_ = column_model.step(one, 0.01, np.abs(one), 1e-9, 1e-4)
        box_state, *_ = box_model.step(box_state, 0.01, np.abs(box_state), 1e-9, 1e-4)
    assert np.all(box_state == box_state[:, :1, :1, :])  # no lateral flux arises from nothing
    peak = np.maximum(np.abs(one).max(axis=1), 1e-12)[:, None]
    # Both systems are small enough for multigrid to solve exactly: rounding remains.
    assert np.max(np.abs(box_state[:, 0, 0, :] - one) / peak) < 1e-10


def test_the_ledger_balances_what_crossed_the_top_to_rounding():
    result = run_in_space(experiment_from_dict(spatial(voxels=(4, 4, 12), duration_h=0.1)))
    balance = result.manifest.outputs["balance"]
    for quantity in ("carbon", "nitrogen", "electrons"):
        assert balance[quantity]["largest_relative_residual"] < 1e-12, quantity
        assert balance[quantity]["imported"] != 0.0  # the box is open, and it shows
    entered = result.manifest.outputs["imported_mol_per_m2"]
    assert entered["oxygen"] > 0
    assert entered["carbon_dioxide"] < 0


def test_a_planted_leak_is_caught_at_the_first_step(monkeypatch):
    honest = ReactionTransport._flows

    def leaky(self, *args):
        into, out = honest(self, *args)
        return into * (1.0 - 1e-3), out  # a thousandth of every arrival goes missing

    monkeypatch.setattr(ReactionTransport, "_flows", leaky)
    with pytest.raises(ConservationError, match="at step 1"):
        run_in_space(experiment_from_dict(spatial(voxels=(12,), duration_h=0.1)))


def test_hour_long_steps_over_an_anoxic_base_stay_positive_and_conserve():
    # 100 um of dense biomass: oxygen penetrates about 60 um, so the base is anoxic.
    thick = [
        {
            "component": "heterotroph",
            "center_um": [],
            "radius_um": 100.0,
            "concentration_mol_per_m3": 2000.0,
        }
    ]
    config = experiment_from_dict(spatial(voxels=(32,), duration_h=1.0, colonies=thick))
    model = build_model(config)
    names = config.network.component_names
    oxygen = names.index("oxygen")
    state = config.domain.initial_state(names, config.initial_mol_per_m3, 0)
    # Nothing is clipped: a clipped value would unbalance the ledger.
    ledger = Ledger(config.network.composition_matrix(), state, QUANTITIES, tolerance=1e-12)
    state, imports, _ = model.integrate(
        state, 0.05, relative_tolerance=1e-3, absolute_tolerance=1e-9
    )
    ledger.exchange(imports)
    ledger.check(state, step=0, time_h=0.05)
    assert state[oxygen, 0] < 1e-3  # the base is anoxic before the long steps begin
    for step in range(1, 4):  # hour-long: far beyond every diffusion time
        state, imports, *_ = model.step(state, 1.0, np.abs(state), 1e-9, 1e-4)
        assert np.all(state >= 0)
        ledger.exchange(imports)
        ledger.check(state, step=step, time_h=0.05 + step)
    assert state[oxygen, 0] < 1e-3


@pytest.mark.parametrize(("threshold", "rounds"), [(None, 1), (0.0, 51)])
def test_the_limiter_keeps_traces_below_the_smallest_normal_number_positive(
    monkeypatch, threshold, rounds
):
    # A voxel holding 4 units in the last place sends one through each of its six
    # faces. Scaled by 4/6, each transfer rounds back up to one unit, so no relative
    # scaling brings the outflow down to the holding: the 64 x 64 x 32 benchmark
    # failed on such traces. With the threshold at zero, the rounds stall and the
    # fallback's last pass must settle it instead.
    import marse.core.implicit as implicit

    if threshold is not None:
        monkeypatch.setattr(implicit, "_TINY", threshold)
    unit = np.finfo(float).smallest_subnormal
    model = ReactionTransport(
        Diffusion(Grid((3, 3, 3), 1.0), np.array([1.0]), np.array([0.0])),
        np.zeros((0, 1)),
        lambda c: np.zeros((0, *c.shape[1:])),
        lambda c: np.zeros((0, 1, *c.shape[1:])),
    )
    y = np.zeros((1, 3, 3, 3))
    y[0, 1, 1, 1] = 4 * unit
    transfers = tuple(np.zeros((1, 3, 3, 3)) for _ in range(3))
    for axis in range(3):
        below = [0, 1, 1, 1]
        below[axis + 1] = 0
        transfers[axis][0, 1, 1, 1] = unit  # out through the upper face
        transfers[axis][tuple(below)] = -unit  # and through the lower one
    extents = np.zeros((0, 3, 3, 3))
    assert model._apply(y, transfers, extents)[0][0, 1, 1, 1] < 0
    new, *_, used = model._limited_update(y, transfers, extents)
    assert np.all(new >= 0)
    assert new.sum() == y.sum()  # nothing clipped
    assert used == rounds


def test_the_run_replays_bit_for_bit_and_an_edit_is_refused(tmp_path):
    config = experiment_from_dict(spatial(voxels=(4, 12), duration_h=0.1))
    first = run_in_space(config)
    second = run_in_space(config)
    np.testing.assert_array_equal(first.final_state, second.final_state)
    path = first.manifest.write(tmp_path / "manifest.json")
    manifest = Manifest.read(path)
    assert manifest.kind == "reactive_transport"
    assert manifest.experiment() == config
    edited = json.loads(path.read_text("utf-8"))
    edited["config"]["domain"]["voxel_um"] = 5.0
    path.write_text(json.dumps(edited), "utf-8")
    with pytest.raises(ValueError, match="does not match its checksum"):
        Manifest.read(path).experiment()


@pytest.mark.invariance
@pytest.mark.parametrize(
    ("scene", "hours"),
    [
        ("networks/surface_biofilm_1d.json", "1.0"),
        ("networks/spreading_column.json", "1.0"),  # biomass spread every quarter hour
        ("environments/oral/stephan_rinse.json", "0.05"),  # the rinse taken, held and spat out
    ],
)
def test_a_run_gives_the_same_result_on_any_number_of_threads(scene, hours):
    # OpenBLAS factorises a matrix of more than about 100 unknowns in parallel, and
    # its last bits then depend on the thread count. The column's coarsest level has
    # 112, so before the coarse solve left LAPACK this run replayed only on a machine
    # with the same number of threads.
    script = (
        "import json, sys\n"
        "from marse.core.reactive_transport import run\n"
        "from marse.schemas.experiment import experiment_from_dict\n"
        "raw = json.loads(open(sys.argv[1], encoding='utf-8').read())\n"
        "raw['duration_h'] = float(sys.argv[2])\n"
        "print(run(experiment_from_dict(raw)).manifest.outputs['final_state_sha256'])\n"
    )
    column = ROOT / "examples" / scene
    digests = set()
    for threads in ("1", "4"):
        variables = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
        env = dict(os.environ, **dict.fromkeys(variables, threads))
        result = subprocess.run(
            [sys.executable, "-c", script, str(column), hours],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        digests.add(result.stdout.strip())
    assert len(digests) == 1, digests


def test_the_rate_jacobian_matches_finite_differences_of_what_the_engine_evaluates():
    config = experiment_from_dict(spatial())
    terms = compile_rates(config.network)
    rng = np.random.default_rng(0)
    scale = np.array([5, 0.2, 5, 1, 3, 800, 400], dtype=float)[:, None, None]
    c = rng.uniform(0.001, 2.0, size=(7, 5, 4)) * scale
    c[rng.random(c.shape) < 0.3] *= -1  # stage values that undershoot: the rates see zero

    def rates(x):
        return process_rates(terms, np.maximum(x, 0.0))

    jacobian = rate_jacobian(terms, c)
    for k in range(7):
        d = np.zeros_like(c)
        d[k] = 1e-6 * np.maximum(1.0, np.abs(c[k]))
        numerical = (rates(c + d) - rates(c - d)) / (2 * d[k])
        np.testing.assert_allclose(
            jacobian[:, k], numerical, rtol=1e-6, atol=1e-9 * np.abs(jacobian).max()
        )
    assert np.all(jacobian[:, c < 0] == 0)  # no slope below zero, so Newton is not misled


# --- the command line -----------------------------------------------------------------


def test_marse_leaves_numpy_error_handling_as_it_found_it(capsys):
    # It used to set every floating-point error to "warn" on the way out, so one
    # test that ran the command line turned underflow into an error for the rest.
    with np.errstate(under="raise"):
        assert main(["check", str(NETWORK)]) == 0
        assert np.geterr()["under"] == "raise"
    assert "glucose_cross_feeding.json" in capsys.readouterr().out


def test_marse_run_writes_totals_frames_and_paraview_files_and_replays(tmp_path, capsys):
    path = tmp_path / "surface.json"
    path.write_text(json.dumps(spatial(voxels=(4, 4, 12), duration_h=0.1)), "utf-8")
    output = tmp_path / "run"
    assert main(["run", str(path), "-o", str(output)]) == 0
    out = capsys.readouterr().out
    assert "kind        reactive_transport" in out
    assert "3-D, 4 x 4 x 12 voxels of 4 um" in out
    assert (output / "totals.csv").read_text("utf-8").startswith("time_h,glucose_mol_per_m2,")
    assert (output / "frames" / "index.json").is_file()
    assert (output / "vtk" / "run.pvd").is_file()
    assert (output / "vtk" / "frame_0000.vti").is_file()
    assert main(["replay", str(output / "manifest.json")]) == 0
    assert "reproduced the recorded results exactly" in capsys.readouterr().out


def test_marse_check_describes_the_space_a_network_will_run_in(tmp_path, capsys):
    path = tmp_path / "surface.json"
    path.write_text(json.dumps(spatial(voxels=(8, 8, 16))), "utf-8")
    assert main(["check", str(path)]) == 0
    out = capsys.readouterr().out
    assert "3-D, 8 x 8 x 16 voxels of 4 um (32 x 32 x 64 um)" in out
    assert "an explicit step would have to be at most" in out
    assert "runnable" in out
