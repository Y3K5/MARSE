"""The well-mixed engine: a reaction network running in a closed box.

The engine is checked against the analytical Monod batch solution (case V2).
It must hold carbon, nitrogen and electrons over ten thousand steps. The
requirements behind known defects 2, 4, 5 and 7 are tested on it directly,
and it must replay bit for bit from its manifest.
"""

import csv
import json
import math
from pathlib import Path

import numpy as np
import pytest

from marse.cli import main
from marse.core.integrators import positive_conservative_step
from marse.core.provenance import Manifest
from marse.core.simulation import ConservationError
from marse.core.well_mixed import ENGINE_VERSION, ExhaustedComponentError, run
from marse.microbes.kinetics import compile_rates, process_rates
from marse.schemas import experiment_from_dict
from marse.validation.analytical import batch_final_biomass, monod_batch_time

EXAMPLE = (
    Path(__file__).resolve().parents[1] / "examples" / "networks" / "glucose_cross_feeding.json"
)

MU, K, Y, S0, X0 = 0.5, 0.3, 1.5, 2.0, 0.05


def yeast(duration_h=6.0, timestep_h=0.5, **run_settings) -> dict:
    """Monod growth on glucose alone: glucose is the only thing consumed.

    The biomass holds no nitrogen (CH1.8O0.5), so the balances close with
    ethanol and carbon dioxide, both produced. The system is then exactly the
    analytical Monod batch problem, with X + Y S conserved.
    """
    return {
        "schema_version": 2,
        "experiment_id": "yeast",
        "components": [
            {"name": "glucose", "phase": "dissolved", "formula": "C6H12O6"},
            {"name": "ethanol", "phase": "dissolved", "formula": "C2H6O"},
            {"name": "carbon_dioxide", "phase": "dissolved", "formula": "CO2"},
            {"name": "yeast", "phase": "particulate", "formula": "CH1.8O0.5"},
        ],
        "processes": [
            {
                "name": "growth",
                "kind": "growth",
                "biomass": "yeast",
                "substrate": "glucose",
                "yield_mol_per_mol": Y,
                "balanced_by": ["ethanol", "carbon_dioxide"],
                "rate": {
                    "maximum_per_h": MU,
                    "factors": [
                        {"component": "glucose", "form": "monod", "half_saturation_mol_per_m3": K}
                    ],
                },
            }
        ],
        "initial_mol_per_m3": {"glucose": S0, "yeast": X0},
        "duration_h": duration_h,
        "timestep_h": timestep_h,
        **run_settings,
    }


def glucose_at(t: float) -> float:
    """Invert the analytical t(S) by bisection: S at time t."""
    low, high = 1e-12, S0
    for _ in range(200):
        middle = 0.5 * (low + high)
        if monod_batch_time(middle, S0, X0, MU, K, Y) > t:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def final(result) -> dict[str, float]:
    return result.final_mol_per_m3


# --- against analytical solutions (case V2) -------------------------------


@pytest.mark.numerical
def test_monod_batch_growth_matches_the_analytical_solution():
    result = run(experiment_from_dict(yeast(relative_tolerance=1e-9)))
    exact = glucose_at(6.0)
    assert final(result)["glucose"] == pytest.approx(exact, rel=1e-6)
    assert final(result)["yeast"] == pytest.approx(X0 + Y * (S0 - exact), rel=1e-6)


@pytest.mark.numerical
def test_the_integration_is_second_order_on_monod_batch_growth():
    config = experiment_from_dict(yeast())
    stoichiometry = config.network.stoichiometric_matrix()
    terms = compile_rates(config.network)
    exact = glucose_at(6.0)
    errors = []
    for n in (100, 200, 400, 800):
        y = np.array([S0, 0.0, 0.0, X0])
        for _ in range(n):
            y, _ = positive_conservative_step(
                y, 6.0 / n, lambda c: process_rates(terms, c), stoichiometry
            )
        errors.append(abs(y[0] - exact))
    orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert np.all(orders > 1.9), orders


def test_the_culture_ends_at_the_biomass_the_substrate_can_pay_for():
    result = run(experiment_from_dict(yeast(duration_h=60.0, timestep_h=1.0)))
    assert final(result)["yeast"] == pytest.approx(batch_final_biomass(S0, X0, Y), rel=1e-6)
    assert final(result)["glucose"] < 1e-6 * S0


# --- the balance ------------------------------------------------------------------


def test_a_closed_box_holds_carbon_nitrogen_and_electrons_over_ten_thousand_steps():
    raw = json.loads(EXAMPLE.read_text("utf-8"))
    raw |= {"duration_h": 100.0, "timestep_h": 0.01, "record_interval_h": 10.0}
    config = experiment_from_dict(raw)
    assert config.steps == 10_000
    balance = run(config).manifest.outputs["balance"]
    for quantity in ("carbon", "nitrogen", "electrons"):
        assert balance[quantity]["largest_relative_residual"] <= 1e-12, quantity


def test_a_leak_is_caught_and_stops_the_run(monkeypatch):
    # The ledger has to be able to fail: plant a leak and it must refuse the run.
    import marse.core.well_mixed as engine

    real = engine.integrate

    def leaking(state, *args, **kwargs):
        state, stats = real(state, *args, **kwargs)
        state = state.copy()
        state[0] += 1e-6  # glucose from nowhere
        return state, stats

    monkeypatch.setattr(engine, "integrate", leaking)
    with pytest.raises(ConservationError, match="carbon is not conserved at step 1"):
        run(experiment_from_dict(json.loads(EXAMPLE.read_text("utf-8"))))


# --- the requirements behind the known defects -----------------------------------


def test_nothing_is_consumed_when_nothing_grows():
    # Known defect 2: the v1 engine consumed substrate at potential, not actual, growth.
    raw = yeast()
    raw["processes"][0]["rate"]["maximum_per_h"] = 0.0
    result = run(experiment_from_dict(raw))
    np.testing.assert_array_equal(
        result.concentrations_mol_per_m3[-1], result.concentrations_mol_per_m3[0]
    )


def test_consumption_is_growth_divided_by_the_yield():
    # Known defect 2: consumption must follow growth, at the stated yield, exactly.
    rows = run(experiment_from_dict(yeast())).concentrations_mol_per_m3
    grown = rows[:, 3] - X0
    consumed = S0 - rows[:, 0]
    np.testing.assert_allclose(consumed, grown / Y, rtol=1e-12, atol=1e-15)


def test_a_monod_factor_halves_the_rate_at_its_half_saturation():
    # Known defect 4: the v1 engine applied the substrate term twice (a quarter here).
    config = experiment_from_dict(yeast())
    rates = process_rates(compile_rates(config.network), np.array([K, 0.0, 0.0, X0]))
    assert rates[0] == pytest.approx(MU * X0 / 2, rel=1e-15)


def test_a_step_far_too_large_stays_positive_and_conserves():
    # Known defect 5: the v1 engine clipped negatives, creating matter.
    config = experiment_from_dict(json.loads(EXAMPLE.read_text("utf-8")))
    network = config.network
    terms = compile_rates(network)
    y = np.array([config.initial_mol_per_m3[n] for n in network.component_names])
    composition = network.composition_matrix()
    before = composition.T @ y
    for _ in range(5):
        y, _ = positive_conservative_step(
            y, 1e3, lambda c: process_rates(terms, c), network.stoichiometric_matrix()
        )
        assert np.all(y >= 0)
    np.testing.assert_allclose(composition.T @ y, before, rtol=1e-13)


def reordered(raw: dict) -> dict:
    swapped = json.loads(json.dumps(raw))
    swapped["components"].reverse()
    swapped["processes"].reverse()
    return swapped


@pytest.mark.invariance
def test_the_order_of_components_and_processes_does_not_matter():
    # Known defect 7: the v1 engine updated species one after another, so listing
    # them in another order changed results by 0.9%. Here every process is
    # updated from the same state: fixed steps agree to rounding, and adaptive
    # runs to well within their tolerance.
    raw = json.loads(EXAMPLE.read_text("utf-8"))
    results = []
    for document in (raw, reordered(raw)):
        config = experiment_from_dict(document)
        network = config.network
        terms = compile_rates(network)
        y = np.array([config.initial_mol_per_m3[n] for n in network.component_names])
        for _ in range(40):
            y, _ = positive_conservative_step(
                y, 0.05, lambda c, t=terms: process_rates(t, c), network.stoichiometric_matrix()
            )
        fixed = dict(zip(network.component_names, y, strict=True))
        results.append((fixed, final(run(config))))
    (fixed, adaptive), (fixed_swapped, adaptive_swapped) = results
    for name in fixed:
        assert fixed_swapped[name] == pytest.approx(fixed[name], rel=1e-12, abs=1e-15), name
        assert adaptive_swapped[name] == pytest.approx(adaptive[name], rel=1e-5, abs=1e-9), name


def test_a_component_assumed_in_excess_that_runs_out_stops_the_run():
    raw = json.loads(EXAMPLE.read_text("utf-8"))
    growth = raw["processes"][0]
    growth["rate"]["factors"] = growth["rate"]["factors"][:2]  # drop the ammonium factor
    growth["rate"]["assumed_in_excess"] = ["ammonium"]
    raw["initial_mol_per_m3"]["ammonium"] = 0.001
    with pytest.raises(ExhaustedComponentError, match="'ammonium' ran out in the step from t = "):
        run(experiment_from_dict(raw))


# --- the example, the record and the command line ------------------------------------


def test_the_example_runs_out_of_oxygen_and_turns_to_fermentation():
    result = run(experiment_from_dict(json.loads(EXAMPLE.read_text("utf-8"))))
    end = final(result)
    assert end["oxygen"] < 1e-9
    assert end["glucose"] < 1e-9
    assert end["lactate"] > 30.0  # most of the glucose became lactate
    assert end["fermenter"] > 10 * end["heterotroph"]
    assert all(v >= 0 for v in result.concentrations_mol_per_m3.ravel())


def test_the_manifest_records_what_ran_and_replays_exactly(tmp_path):
    result = run(experiment_from_dict(json.loads(EXAMPLE.read_text("utf-8"))))
    manifest = result.manifest
    assert manifest.kind == "well_mixed"
    assert manifest.models["engine"] == ENGINE_VERSION
    assert manifest.random_streams == ()
    assert set(manifest.outputs) == {
        "final_time_h",
        "final_mol_per_m3",
        "balance",
        "substeps",
        "final_state_sha256",
    }
    replayed = run(Manifest.read(manifest.write(tmp_path / "manifest.json")).experiment())
    assert replayed.manifest.outputs == manifest.outputs
    assert replayed.manifest.run_id == manifest.run_id


def test_an_edited_manifest_is_refused(tmp_path):
    path = run(experiment_from_dict(yeast())).manifest.write(tmp_path / "manifest.json")
    raw = json.loads(path.read_text("utf-8"))
    raw["config"]["initial_mol_per_m3"]["glucose"] = 3.0
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        Manifest.read(path).experiment()


def test_the_trajectory_has_units_in_its_header_and_starts_at_the_initial_state(tmp_path):
    result = run(experiment_from_dict(yeast()))
    path = result.write_trajectory(tmp_path / "trajectory.csv")
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    assert rows[0] == [
        "time_h",
        "glucose_mol_per_m3",
        "ethanol_mol_per_m3",
        "carbon_dioxide_mol_per_m3",
        "yeast_mol_per_m3",
    ]
    assert [float(v) for v in rows[1]] == [0.0, S0, 0.0, 0.0, X0]
    assert len(rows) - 1 == len(result.times_h) == math.ceil(6.0 / 0.5) + 1


def test_marse_run_and_replay_a_network(tmp_path, capsys):
    output = tmp_path / "run"
    assert main(["run", str(EXAMPLE), "-o", str(output)]) == 0
    out = capsys.readouterr().out
    assert "kind        well_mixed" in out
    assert "carbon, nitrogen and electrons conserved to" in out
    assert (output / "trajectory.csv").is_file()
    assert main(["replay", str(output / "manifest.json")]) == 0
    assert "reproduced the recorded results exactly" in capsys.readouterr().out


def test_marse_run_refuses_a_network_that_cannot_run(tmp_path, capsys):
    raw = json.loads(EXAMPLE.read_text("utf-8"))
    del raw["processes"][2]["rate"]
    path = tmp_path / "unrated.json"
    path.write_text(json.dumps(raw), encoding="utf-8")
    assert main(["run", str(path)]) == 2
    assert "every process needs a rate to run" in capsys.readouterr().out


def test_marse_check_says_whether_a_network_can_run(capsys):
    assert main(["check", str(EXAMPLE)]) == 0
    out = capsys.readouterr().out
    assert "rate: 0.4 /h x heterotroph x monod(glucose; K 0.05)" in out
    # The suggested command uses the path as given, so it runs from where check ran.
    assert f"runnable    48 h in steps of at most 0.5 h: marse run {EXAMPLE}" in out
