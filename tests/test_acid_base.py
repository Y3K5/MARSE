"""pH from electroneutrality, and rates that depend on it.

Three of the criteria set before Stage S1 was built are checked here:

- G1: the hydrogen ion concentration agrees with an independent bisection of
  the charge balance to 1e-10;
- G2: every voxel is left neutral to 1e-12 of the charges present;
- G4: the rate Jacobian, through dpH/dc, agrees with finite differences to
  1e-6.
"""

from fractions import Fraction

import numpy as np
import pytest

from marse.chemistry import ChargeBalance, ph_of_hydrogen
from marse.core.config import ConfigError
from marse.microbes.cardinal import cardinal_ph, cardinal_ph_slope
from marse.microbes.kinetics import compile_rates, process_rates, rate_jacobian
from marse.schemas.network import network_from_dict

SALIVA_AND_PLAQUE = {
    "schema_version": 2,
    "pkw": 13.6,
    "components": [
        {"name": "glucose", "phase": "dissolved", "formula": "C6H12O6"},
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
        {
            "name": "ammonium",
            "phase": "dissolved",
            "formula": "NH4",
            "charge": 1,
            "acid_base": {"pka": [9.0]},
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
            "stoichiometry_mol_per_mol": {"glucose": -1, "lactate": 2},
            "rate": {
                "maximum_per_h": 0.5,
                "proportional_to": "bacteria",
                "factors": [
                    {"component": "glucose", "form": "monod", "half_saturation_mol_per_m3": 1},
                    {"form": "ph", "ph_min": 4, "ph_optimum": 7, "ph_max": 9},
                ],
            },
        }
    ],
}
NETWORK = network_from_dict(SALIVA_AND_PLAQUE)
NAMES = NETWORK.component_names
BALANCE = ChargeBalance.of(NETWORK)
TYPICAL = np.array([30, 40, 5, 10, 3, 25, 25, 120, 50, 800], dtype=float)  # plaque, mol/m3


def _column(values):
    return np.asarray(values, dtype=float)[:, None]


def _bisect(balance, c):
    """The root of the charge balance in every voxel, by bisection in log h to the last bit."""
    low = np.full(c.shape[1:], 1e-14)
    high = np.full(c.shape[1:], 1e4)
    for _ in range(400):
        middle = np.sqrt(low * high)
        net, _ = balance.residual(c, middle)
        low = np.where(net < 0, middle, low)
        high = np.where(net < 0, high, middle)
    return np.sqrt(low * high)


# --- G1 and G2 ------------------------------------------------------------------------


@pytest.mark.parametrize("acid", ["lactate", "carbonate", "phosphate", "ammonium"])
def test_g1_the_charge_balance_is_solved_as_bisection_would(acid):
    rng = np.random.default_rng(1)
    c = np.zeros((len(NAMES), 500))
    c[NAMES.index(acid)] = rng.uniform(0.1, 50, 500)
    c[NAMES.index("potassium")] = rng.uniform(0, 60, 500)
    c[NAMES.index("chloride")] = rng.uniform(0, 60, 500)
    h = BALANCE.hydrogen(c)
    assert np.max(np.abs(h - _bisect(BALANCE, c)) / h) < 1e-10


def test_g2_every_voxel_is_left_neutral():
    rng = np.random.default_rng(2)
    c = np.abs(rng.normal(size=(len(NAMES), 2000))) * TYPICAL[:, None]
    h = BALANCE.hydrogen(c)
    net, _ = BALANCE.residual(c, h)
    ions = np.abs(c[[NAMES.index("potassium"), NAMES.index("chloride")]]).sum(axis=0) + h
    assert np.max(np.abs(net) / ions) < 1e-12
    assert ph_of_hydrogen(h).min() > 1.0  # the range is exercised
    assert ph_of_hydrogen(h).max() < 13.0


def test_a_weak_acid_with_its_salt_gives_back_the_ph_it_was_made_at():
    """Inverse check: the potassium that neutralises lactate at a chosen pH gives back that pH."""
    kw = 10**-13.6 * 1e6
    for ph in (3.0, 3.86, 5.0, 7.0, 8.5):
        h = 1000.0 * 10**-ph
        total = 20.0
        ka = 1000.0 * 10**-3.86
        dissociated = total * ka / (ka + h)
        c = np.zeros((len(NAMES), 1))
        c[NAMES.index("lactate")] = total
        c[NAMES.index("potassium")] = dissociated - h + kw / h
        assert ph_of_hydrogen(BALANCE.hydrogen(c))[0] == pytest.approx(ph, abs=1e-12)


def test_at_the_pka_a_buffer_is_half_dissociated():
    c = np.zeros((len(NAMES), 1))
    c[NAMES.index("lactate")] = 100.0
    c[NAMES.index("potassium")] = 50.0  # half neutralised: Henderson-Hasselbalch gives pKa
    assert ph_of_hydrogen(BALANCE.hydrogen(c))[0] == pytest.approx(3.86, abs=0.01)


def test_pure_water_is_neutral_at_half_pkw():
    c = np.zeros((len(NAMES), 3))
    assert ph_of_hydrogen(BALANCE.hydrogen(c)) == pytest.approx(6.8, abs=1e-12)


def test_the_solve_works_for_a_well_mixed_box_too():
    h = BALANCE.hydrogen(TYPICAL)
    assert h.shape == ()
    np.testing.assert_allclose(h, BALANCE.hydrogen(TYPICAL[:, None])[0], rtol=1e-15)


def test_a_negative_concentration_counts_as_zero_and_moves_nothing():
    c = TYPICAL[:, None].copy()
    c[NAMES.index("chloride")] = -1e-6
    zeroed = c.copy()
    zeroed[NAMES.index("chloride")] = 0.0
    h = BALANCE.hydrogen(c)
    assert h == BALANCE.hydrogen(zeroed)
    assert BALANCE.hydrogen_slopes(c, h)[NAMES.index("chloride")] == 0.0


# --- the slopes and G4 ----------------------------------------------------------------


def _central(k, c, step):
    up, down = c.copy(), c.copy()
    up[k] += step
    down[k] -= step
    return (BALANCE.hydrogen(up) - BALANCE.hydrogen(down)) / (2 * step)


def test_the_slopes_of_h_are_those_of_finite_differences():
    """Richardson-extrapolated central differences, which are fourth order.

    A solved h carries rounding of about 1e-15 h, so plain central differences
    cannot check a slope to better than about 1e-6; extrapolation removes their
    step-squared error and allows a larger step.
    """
    rng = np.random.default_rng(3)
    c = np.abs(rng.normal(size=(len(NAMES), 50))) * TYPICAL[:, None] + 0.1
    h = BALANCE.hydrogen(c)
    slopes = BALANCE.hydrogen_slopes(c, h)
    for k in range(len(NAMES)):
        step = 1e-3 * np.maximum(1.0, c[k])
        numeric = (4 * _central(k, c, step / 2) - _central(k, c, step)) / 3
        noise = 1e-14 * h / step
        np.testing.assert_allclose(slopes[k], numeric, rtol=1e-7, atol=noise.max())


def test_acid_raises_h_and_base_lowers_it():
    h = BALANCE.hydrogen(TYPICAL[:, None])
    slopes = BALANCE.hydrogen_slopes(TYPICAL[:, None], h)
    assert slopes[NAMES.index("lactate")] > 0
    assert slopes[NAMES.index("chloride")] > 0
    assert slopes[NAMES.index("potassium")] < 0
    assert slopes[NAMES.index("glucose")] == 0  # uncharged at any pH


def test_the_cardinal_ph_slope_is_the_derivative_of_the_factor():
    ph = np.linspace(4.05, 8.95, 200)
    step = 1e-6
    numeric = (cardinal_ph(ph + step, 4, 7, 9) - cardinal_ph(ph - step, 4, 7, 9)) / (2 * step)
    np.testing.assert_allclose(cardinal_ph_slope(ph, 4, 7, 9), numeric, rtol=1e-6, atol=1e-9)
    assert cardinal_ph_slope(7.0, 4, 7, 9) == 0.0
    assert cardinal_ph_slope(np.array([3.0, 4.0, 9.0, 10.0]), 4, 7, 9).tolist() == [0, 0, 0, 0]


def test_g4_the_rate_jacobian_follows_the_ph_through_every_charged_component():
    terms = compile_rates(NETWORK)
    rng = np.random.default_rng(4)
    c = np.abs(rng.normal(size=(len(NAMES), 50))) * TYPICAL[:, None] + 0.1
    jacobian = rate_jacobian(terms, c)
    numeric = np.zeros_like(jacobian)
    for k in range(len(NAMES)):
        step = 1e-5 * np.maximum(1.0, c[k])
        up, down = c.copy(), c.copy()
        up[k] += step
        down[k] -= step
        numeric[:, k] = (process_rates(terms, up) - process_rates(terms, down)) / (2 * step)
    scale = np.abs(numeric).max()
    assert np.max(np.abs(jacobian - numeric) / (np.abs(numeric) + 1e-9 * scale)) < 1e-6
    assert np.any(jacobian[0, NAMES.index("lactate")] < 0)  # acid slows its own production


def test_a_ph_factor_scales_the_rate_by_the_cardinal_model():
    terms = compile_rates(NETWORK)
    c = TYPICAL[:, None]
    ph = ph_of_hydrogen(BALANCE.hydrogen(c))
    monod = c[0] / (1.0 + c[0])
    expected = 0.5 * c[NAMES.index("bacteria")] * monod * cardinal_ph(ph, 4, 7, 9)
    np.testing.assert_allclose(process_rates(terms, c)[0], expected, rtol=1e-14)


# --- reading acids, bases and pH factors ----------------------------------------------


def _with(change):
    raw = {
        **SALIVA_AND_PLAQUE,
        "components": [dict(c) for c in SALIVA_AND_PLAQUE["components"]],
        "processes": [dict(p) for p in SALIVA_AND_PLAQUE["processes"]],
    }
    change(raw)
    return raw


def _set_pka(values):
    def change(raw):
        raw["components"][1]["acid_base"] = {"pka": values}

    return change


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (_set_pka([]), "list at least one pKa"),
        (_set_pka([6.8, 2.0]), "must rise strictly"),
        (_set_pka([3.86, 5, 7, 9, 11, 12, 13]), "cannot leave C3H6O3"),
        (_set_pka([20]), "between -2 and 16"),
        (lambda raw: raw.update(pkw=20), "between 11 and 16"),
    ],
)
def test_acid_base_components_are_checked(change, message):
    with pytest.raises(ConfigError, match=message):
        network_from_dict(_with(change))


def _factor(**fields):
    def change(raw):
        rate = dict(raw["processes"][0]["rate"])
        rate["factors"] = [rate["factors"][0], fields]
        raw["processes"][0] = {**raw["processes"][0], "rate": rate}

    return change


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (_factor(form="ph", ph_min=4, ph_optimum=7), "needs ph_max"),
        (_factor(form="ph", ph_min=7, ph_optimum=4, ph_max=9), "ph_min < ph_optimum < ph_max"),
        (
            _factor(form="ph", component="lactate", ph_min=4, ph_optimum=7, ph_max=9),
            "remove component",
        ),
        (
            _factor(form="ph", ph_min=4, ph_optimum=7, ph_max=9, half_saturation_mol_per_m3=1),
            "does not apply to a ph factor",
        ),
        (_factor(form="monod", half_saturation_mol_per_m3=1), "needs component"),
        (
            _factor(form="monod", component="lactate", half_saturation_mol_per_m3=1, ph_min=4),
            "applies only to a ph factor",
        ),
    ],
)
def test_ph_factors_are_checked(change, message):
    with pytest.raises(ConfigError, match=message):
        network_from_dict(_with(change))


def test_a_rate_takes_one_ph_factor():
    def twice(raw):
        rate = dict(raw["processes"][0]["rate"])
        rate["factors"] = [*rate["factors"], rate["factors"][1]]
        raw["processes"][0] = {**raw["processes"][0], "rate": rate}

    with pytest.raises(ConfigError, match="at most one ph factor"):
        network_from_dict(_with(twice))


def test_a_ph_needs_charges_and_pkw_needs_a_ph():
    neutral = {
        "schema_version": 2,
        "components": [
            {"name": "glucose", "phase": "dissolved", "formula": "C6H12O6"},
            {"name": "bacteria", "phase": "particulate", "formula": "CH1.8O0.5N0.2"},
        ],
        "processes": [],
    }
    with pytest.raises(ConfigError, match="pkw: applies only"):
        network_from_dict({**neutral, "pkw": 13.6})
    lactic_acid = {"name": "lactic_acid", "phase": "dissolved", "formula": "C3H6O3"}
    fermenting = {
        **neutral,
        "components": [*neutral["components"], lactic_acid],
        "processes": [
            {
                "name": "fermentation",
                "kind": "reaction",
                "stoichiometry_mol_per_mol": {"glucose": -1, "lactic_acid": 2},
                "rate": {
                    "maximum_per_h": 1,
                    "proportional_to": "bacteria",
                    "factors": [
                        {"component": "glucose", "form": "monod", "half_saturation_mol_per_m3": 1},
                        {"form": "ph", "ph_min": 4, "ph_optimum": 7, "ph_max": 9},
                    ],
                },
            }
        ],
    }
    with pytest.raises(ConfigError, match="needs charged components"):
        network_from_dict(fermenting)


def test_acids_bases_and_ph_factors_are_written_back_as_read():
    written = NETWORK.to_dict()
    assert written["pkw"] == 13.6
    assert NETWORK.pkw == Fraction(68, 5)
    assert written["components"][2]["acid_base"] == {"pka": [6.1, 10]}
    assert written["processes"][0]["rate"]["factors"][1] == {
        "form": "ph",
        "ph_min": 4,
        "ph_optimum": 7,
        "ph_max": 9,
    }
    assert network_from_dict(written) == NETWORK
    assert NETWORK.has_ph
    assert NETWORK.process("fermentation").rate.factors[1].describe() == "ph(4, 7, 9)"


# --- runs that set a pH -------------------------------------------------------------------

import copy  # noqa: E402
import csv  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402

from marse.cli import main  # noqa: E402
from marse.core.reactive_transport import run as run_in_space  # noqa: E402
from marse.core.well_mixed import run as run_well_mixed  # noqa: E402
from marse.schemas import experiment_from_dict  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def _neutral(concentrations, ph=7.0):
    """The potassium that makes a composition neutral at ``ph``, added to what it holds."""
    c = np.array([concentrations.get(n, 0.0) for n in NAMES])
    net, _ = BALANCE.residual(c, np.array(1000.0 * 10**-ph))
    return {**concentrations, "potassium": concentrations.get("potassium", 0.0) - float(net)}


SALIVA = _neutral({"glucose": 50.0, "carbonate": 5.0, "phosphate": 10.0, "chloride": 20.0})
_KA, _H7 = 1000.0 * 10**-4.8, 1000.0 * 10**-7.0
FIXED = {  # carboxyl groups on the bacteria, with the cations bound to them at pH 7
    "carboxyl_groups": 120.0,
    "bound_potassium": 120.0 * _KA / (_KA + _H7),
    "bacteria": 800.0,
}
PLAQUE = {**SALIVA, **FIXED}


def _box():
    raw = copy.deepcopy(SALIVA_AND_PLAQUE)
    raw.update(experiment_id="acid_box", duration_h=1.0, timestep_h=0.1, initial_mol_per_m3=PLAQUE)
    return raw


def test_a_fermenting_box_acidifies_and_records_its_ph(tmp_path):
    result = run_well_mixed(experiment_from_dict(_box()))
    assert result.ph is not None
    assert result.ph[0] == pytest.approx(7.0, abs=1e-9)
    assert np.all(np.diff(result.ph) < 0)
    outputs = result.manifest.outputs
    assert outputs["ph"]["final"] == pytest.approx(result.ph[-1])
    assert outputs["ph"]["lowest_at_h"] == 1.0
    assert list(outputs["balance"]) == [
        "carbon",
        "nitrogen",
        "electrons",
        "phosphorus",
        "potassium",
        "chlorine",
    ]
    with result.write_trajectory(tmp_path / "trajectory.csv").open(encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    assert rows[0][-1] == "ph"
    assert float(rows[-1][-1]) == pytest.approx(result.ph[-1], abs=1e-6)


def _column():
    raw = _box()
    raw.update(experiment_id="acid_column", duration_h=0.25, timestep_h=0.05)
    saliva = SALIVA
    dissolved = [c["name"] for c in raw["components"] if c["phase"] == "dissolved"]
    raw["domain"] = {
        "voxels": [16],
        "voxel_um": 10,
        "bulk_mol_per_m3": saliva,
        "diffusivity_m2_per_s": {n: 1e-9 for n in dissolved},
        "colonies": [
            {
                "component": n,
                "center_um": [],
                "radius_um": 80,
                "concentration_mol_per_m3": PLAQUE[n],
            }
            for n in ("carboxyl_groups", "bound_potassium", "bacteria")
        ],
    }
    raw["initial_mol_per_m3"] = saliva
    return raw


def test_a_column_records_the_ph_at_the_substratum(tmp_path):
    path = tmp_path / "column.json"
    path.write_text(json.dumps(_column()), encoding="utf-8")
    assert main(["run", str(path), "-o", str(tmp_path / "out")]) == 0
    with (tmp_path / "out" / "ph.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert list(rows[0]) == [
        "time_h",
        "substratum_mean",
        "substratum_min",
        "substratum_max",
        "box_min",
        "box_max",
    ]
    assert float(rows[0]["substratum_mean"]) == pytest.approx(7.0, abs=1e-6)
    assert float(rows[-1]["substratum_mean"]) < 6.9  # the plaque at the bottom acidifies
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["outputs"]["ph"]["lowest_at_substratum"] == pytest.approx(
        float(rows[-1]["substratum_min"]), abs=1e-6
    )


def test_marse_check_prints_the_ph_each_composition_implies(tmp_path, capsys):
    path = tmp_path / "column.json"
    path.write_text(json.dumps(_column()), encoding="utf-8")
    assert main(["check", str(path)]) == 0
    out = " ".join(capsys.readouterr().out.split())  # wrapped lines joined
    assert "pH from the charge balance, pKw 13.6" in out
    assert "lactate (3.86)" in out
    assert "phosphate (2 6.8 11.7)" in out
    assert "pH 7.00 in the bulk liquid; 7.00 to 7.00 in the box at the start" in out
    assert "each conserving carbon, nitrogen, electrons, phosphorus, potassium and chlorine" in out


# --- an absolute tolerance per component ------------------------------------------------


def test_a_tolerance_per_component_is_read_completed_and_written_back():
    raw = _box()
    raw["absolute_tolerance_mol_per_m3"] = {"glucose": 1e-6}
    config = experiment_from_dict(raw)
    tolerance = config.absolute_tolerance_mol_per_m3
    assert tolerance["glucose"] == 1e-6
    assert tolerance["lactate"] == 1e-9
    assert list(tolerance) == list(NAMES)
    assert config.absolute_tolerances().tolist() == [tolerance[n] for n in NAMES]
    assert experiment_from_dict(json.loads(json.dumps(config.to_dict()))) == config


@pytest.mark.parametrize(
    ("tolerance", "message"),
    [({"glucos": 1e-6}, "did you mean 'glucose'"), ({"glucose": 0}, "glucose: must be positive")],
)
def test_a_tolerance_per_component_is_checked(tolerance, message):
    raw = _box()
    raw["absolute_tolerance_mol_per_m3"] = tolerance
    with pytest.raises(ConfigError, match=message):
        experiment_from_dict(raw)


def test_one_tolerance_given_per_component_changes_nothing():
    """The same value for every component takes the per-component path and gives the same run."""
    one = _box()
    each = _box()
    each["absolute_tolerance_mol_per_m3"] = {n: 1e-9 for n in NAMES}
    assert (
        run_well_mixed(experiment_from_dict(one)).manifest.outputs["final_state_sha256"]
        == run_well_mixed(experiment_from_dict(each)).manifest.outputs["final_state_sha256"]
    )
    one, each = _column(), _column()
    each["absolute_tolerance_mol_per_m3"] = {n: 1e-9 for n in NAMES}
    assert (
        run_in_space(experiment_from_dict(one)).manifest.outputs["final_state_sha256"]
        == run_in_space(experiment_from_dict(each)).manifest.outputs["final_state_sha256"]
    )
