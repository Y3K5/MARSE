"""Tests for the Phase 1 simulation kernel: config, state, seeds, manifest, run loop."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from marse.core.config import ConfigError, ExperimentConfig, experiment_from_dict, load_experiment
from marse.core.provenance import Manifest, config_checksum
from marse.core.seeds import SeedRegistry
from marse.core.simulation import ConservationError, UnstableStepError, run
from marse.validation.analytical import batch_final_biomass, monod_batch_time

EXAMPLE = (
    Path(__file__).resolve().parents[1] / "examples" / "experiments" / "two_species_batch.json"
)


def minimal_experiment(**overrides) -> dict:
    """A one-organism experiment with no maintenance, decay or environmental scaling."""
    config = {
        "experiment_id": "unit-test",
        "seed": 7,
        "duration_h": 12.0,
        "timestep_h": 0.002,
        "checkpoint_interval_h": 1.0,
        "environment": {"temperature_c": 37.0, "ph": 7.0},
        "substrate": {"name": "glucose", "initial_mm": 10.0},
        "organisms": [
            {
                "name": "solo",
                "initial_biomass_g_per_l": 0.01,
                "mu_opt_per_h": 0.8,
                "k_s_mm": 0.4,
                "yield_g_per_mmol": 0.08,
            }
        ],
    }
    return config | overrides


# --- configuration validation ------------------------------------------------


def test_example_experiment_loads_and_is_valid():
    config = load_experiment(EXAMPLE)
    assert config.experiment_id == "two-species-batch"
    assert [o.name for o in config.organisms] == ["fast_grower", "gleaner"]
    assert config.steps == 24_000


@pytest.mark.parametrize(
    ("mutate", "expected"),
    [
        (lambda c: c.__setitem__("duration_h", 0.0), "duration_h"),
        (lambda c: c.__setitem__("timestep_h", -1.0), "timestep_h"),
        (lambda c: c.__setitem__("seed", -3), "seed"),
        (lambda c: c.__setitem__("organisms", []), "organisms"),
        (lambda c: c.__setitem__("checkpoint_interval_h", 1e-9), "checkpoint_interval_h"),
        (lambda c: c["organisms"][0].__setitem__("k_s_mm", 0.0), "k_s_mm"),
        (lambda c: c["organisms"][0].__setitem__("yield_g_per_mmol", -1.0), "yield_g_per_mmol"),
        (
            lambda c: c["organisms"][0].__setitem__("initial_biomass_g_per_l", 0.0),
            "initial_biomass",
        ),
        (lambda c: c["organisms"][0].__setitem__("decay_per_h", -0.1), "decay_per_h"),
        (lambda c: c["environment"].__setitem__("ph", 20.0), "ph"),
        (lambda c: c["substrate"].__setitem__("initial_mm", -1.0), "initial_mm"),
    ],
)
def test_invalid_configurations_are_rejected_with_a_named_field(mutate, expected):
    config = minimal_experiment()
    mutate(config)
    with pytest.raises(ConfigError, match=expected):
        experiment_from_dict(config)


def test_unknown_fields_are_rejected_rather_than_silently_ignored():
    with pytest.raises(ConfigError, match="unknown field"):
        experiment_from_dict(minimal_experiment(tempreature_c=37.0))
    config = minimal_experiment()
    config["organisms"][0]["mu_max_per_h"] = 1.0  # a plausible-looking typo
    with pytest.raises(ConfigError, match="unknown field"):
        experiment_from_dict(config)


def test_missing_required_fields_are_named():
    config = minimal_experiment()
    del config["substrate"]
    with pytest.raises(ConfigError, match="substrate"):
        experiment_from_dict(config)


def test_duplicate_organism_names_are_rejected():
    config = minimal_experiment()
    config["organisms"].append(dict(config["organisms"][0]))
    with pytest.raises(ConfigError, match="unique"):
        experiment_from_dict(config)


def test_cardinal_temperatures_outside_the_ctmi_domain_are_rejected():
    config = minimal_experiment()
    config["organisms"][0]["cardinal_temperature_c"] = [5.0, 20.0, 45.0]  # opt below the midpoint
    with pytest.raises(ConfigError, match="CTMI"):
        experiment_from_dict(config)


def test_invalid_json_is_reported_with_the_line(tmp_path):
    path = tmp_path / "broken.json"
    path.write_text('{"experiment_id": "x",\n  "seed": }\n')
    with pytest.raises(ConfigError, match="invalid JSON"):
        load_experiment(path)


# --- seeded streams ----------------------------------------------------------


def test_streams_are_reproducible_for_the_same_seed():
    a = SeedRegistry(42).stream("growth").random(5)
    b = SeedRegistry(42).stream("growth").random(5)
    np.testing.assert_array_equal(a, b)


def test_streams_differ_between_names_and_between_seeds():
    registry = SeedRegistry(42)
    assert not np.array_equal(
        registry.stream("growth").random(5), registry.stream("decay").random(5)
    )
    assert not np.array_equal(
        SeedRegistry(1).stream("growth").random(5), SeedRegistry(2).stream("growth").random(5)
    )


def test_a_streams_values_do_not_depend_on_what_else_was_requested():
    """Adding a provider must not perturb the stream any other provider sees."""
    alone = SeedRegistry(99).stream("growth").random(4)
    crowded = SeedRegistry(99)
    crowded.stream("diffusion")
    crowded.stream("detachment")
    np.testing.assert_array_equal(crowded.stream("growth").random(4), alone)


def test_the_same_name_returns_the_same_stream_within_a_run():
    registry = SeedRegistry(5)
    assert registry.stream("growth") is registry.stream("growth")


def test_issued_names_are_sorted_for_a_stable_manifest():
    registry = SeedRegistry(5)
    for name in ("zeta", "alpha", "mu"):
        registry.stream(name)
    assert registry.issued_names == ("alpha", "mu", "zeta")


@pytest.mark.parametrize("bad", [-1, 1.5, True, "7"])
def test_invalid_seeds_are_rejected(bad):
    with pytest.raises(ValueError, match="seed"):
        SeedRegistry(bad)


# --- the run loop ------------------------------------------------------------


def test_single_species_run_matches_the_analytical_batch_solution():
    """The kernel must reproduce the closed-form solution it is built from."""
    config = experiment_from_dict(minimal_experiment())
    result = run(config)
    organism = config.organisms[0]
    mu_max = organism.mu_opt_per_h  # no cardinal scaling in this configuration

    # Compare the time to reach each checkpoint's substrate level with the exact solution.
    checked = 0
    for t, s in zip(result.times_h, result.substrate_mm, strict=True):
        if s <= 1e-3 * config.substrate.initial_mm or t == 0.0:
            continue
        exact = float(
            monod_batch_time(
                s,
                config.substrate.initial_mm,
                organism.initial_biomass_g_per_l,
                mu_max,
                organism.k_s_mm,
                organism.yield_g_per_mmol,
            )
        )
        assert t == pytest.approx(exact, rel=1e-6), f"at S = {s} mM"
        checked += 1
    assert checked >= 3

    expected = batch_final_biomass(
        config.substrate.initial_mm, organism.initial_biomass_g_per_l, organism.yield_g_per_mmol
    )
    assert result.final_state.total_biomass_g_per_l == pytest.approx(expected, rel=1e-6)


def test_substrate_is_conserved_and_never_negative():
    result = run(experiment_from_dict(minimal_experiment()))
    final = result.final_state
    assert final.substrate_mm >= 0.0
    assert final.substrate_mm + final.substrate_consumed_mm == pytest.approx(
        result.config.substrate.initial_mm, abs=1e-9
    )
    assert np.all(result.substrate_mm >= 0.0)
    assert np.all(np.diff(result.substrate_mm) <= 1e-12)  # substrate only ever decreases


def test_a_non_finite_state_stops_the_run(monkeypatch):
    """A blow-up must be reported, not returned as a result."""
    import marse.core.simulation as simulation

    original = simulation._derivatives
    calls = {"n": 0}

    def exploding(*args, **kwargs):
        calls["n"] += 1
        growth, consumption = original(*args, **kwargs)
        if calls["n"] > 20:
            return growth * np.inf, consumption
        return growth, consumption

    monkeypatch.setattr(simulation, "_derivatives", exploding)
    with pytest.raises(ConservationError, match=r"non-finite|diverged"):
        run(experiment_from_dict(minimal_experiment()))


def test_a_timestep_too_large_for_the_rates_is_refused():
    """Silently clamping an unstable run would hide the error; it is raised instead."""
    config = minimal_experiment(duration_h=10.0, timestep_h=5.0, checkpoint_interval_h=5.0)
    config["organisms"][0]["decay_per_h"] = 3.0  # decay per step far exceeds the biomass
    with pytest.raises(UnstableStepError, match="too large"):
        run(experiment_from_dict(config))


def test_a_stable_timestep_is_accepted():
    config = minimal_experiment(duration_h=10.0, timestep_h=0.001)
    config["organisms"][0]["decay_per_h"] = 3.0
    result = run(experiment_from_dict(config))
    assert np.all(result.biomass_g_per_l >= 0.0)


def test_run_ends_exactly_on_the_requested_duration():
    config = experiment_from_dict(minimal_experiment(duration_h=5.0, timestep_h=0.003))
    result = run(config)
    assert result.final_state.time_h == pytest.approx(5.0, abs=1e-9)
    assert result.times_h[-1] == pytest.approx(5.0, abs=1e-9)


def test_decay_reduces_biomass_once_the_substrate_is_gone():
    config = experiment_from_dict(minimal_experiment(duration_h=24.0))
    without = run(config).final_state.total_biomass_g_per_l
    decaying = minimal_experiment(duration_h=24.0)
    decaying["organisms"][0]["decay_per_h"] = 0.02
    assert run(experiment_from_dict(decaying)).final_state.total_biomass_g_per_l < without


def test_temperature_scaling_slows_growth_away_from_the_optimum():
    def final_biomass(temperature: float) -> float:
        config = minimal_experiment(duration_h=4.0)
        config["environment"]["temperature_c"] = temperature
        config["organisms"][0]["cardinal_temperature_c"] = [6.0, 40.0, 47.0]
        return run(experiment_from_dict(config)).final_state.total_biomass_g_per_l

    assert final_biomass(15.0) < final_biomass(30.0) < final_biomass(40.0)


def test_growth_stops_outside_the_cardinal_range():
    config = minimal_experiment(duration_h=4.0)
    config["environment"]["temperature_c"] = 50.0  # above T_max
    config["organisms"][0]["cardinal_temperature_c"] = [6.0, 40.0, 47.0]
    result = run(experiment_from_dict(config))
    assert result.final_state.total_biomass_g_per_l == pytest.approx(0.01)
    assert result.final_state.substrate_mm == pytest.approx(10.0)


def test_two_species_batch_is_won_by_the_faster_grower():
    result = run(load_experiment(EXAMPLE))
    fast, gleaner = result.final_state.biomass_g_per_l
    # In batch culture with abundant substrate, high mu_max wins; the chemostat
    # result of theory.md section 7.2 is the opposite case, not a contradiction.
    assert fast > gleaner


# --- determinism and replay --------------------------------------------------


def test_the_same_configuration_produces_identical_results():
    config = load_experiment(EXAMPLE)
    first, second = run(config), run(config)
    np.testing.assert_array_equal(first.biomass_g_per_l, second.biomass_g_per_l)
    np.testing.assert_array_equal(first.substrate_mm, second.substrate_mm)
    assert first.manifest.run_id == second.manifest.run_id


def test_replaying_from_a_manifest_reproduces_the_run_exactly(tmp_path):
    original = run(load_experiment(EXAMPLE))
    path = original.manifest.write(tmp_path / "manifest.json")

    restored = Manifest.read(path)
    replayed = run(restored.experiment())

    assert replayed.manifest.run_id == original.manifest.run_id
    np.testing.assert_array_equal(replayed.biomass_g_per_l, original.biomass_g_per_l)
    np.testing.assert_array_equal(replayed.substrate_mm, original.substrate_mm)
    assert replayed.manifest.outputs["final_state"] == original.manifest.outputs["final_state"]


def test_an_edited_manifest_is_refused(tmp_path):
    result = run(experiment_from_dict(minimal_experiment()))
    path = result.manifest.write(tmp_path / "manifest.json")
    tampered = json.loads(path.read_text())
    tampered["config"]["substrate"]["initial_mm"] = 999.0
    path.write_text(json.dumps(tampered))

    with pytest.raises(ValueError, match="checksum"):
        Manifest.read(path).experiment()


def test_manifests_from_a_future_schema_are_refused(tmp_path):
    result = run(experiment_from_dict(minimal_experiment()))
    path = result.manifest.write(tmp_path / "manifest.json")
    raw = json.loads(path.read_text())
    raw["manifest_version"] = 99
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="manifest version"):
        Manifest.read(path)


def test_run_id_depends_only_on_the_configuration():
    a = run(experiment_from_dict(minimal_experiment())).manifest.run_id
    b = run(experiment_from_dict(minimal_experiment(seed=8))).manifest.run_id
    assert a != b  # the seed is part of the configuration
    assert a == run(experiment_from_dict(minimal_experiment())).manifest.run_id
    assert a.startswith("MARSE-")


def test_checksum_changes_when_any_parameter_changes():
    base = config_checksum(experiment_from_dict(minimal_experiment()))
    nudged = minimal_experiment()
    nudged["organisms"][0]["k_s_mm"] = 0.4000001
    assert config_checksum(experiment_from_dict(nudged)) != base


# --- the manifest must not identify the person or machine --------------------


def test_manifest_records_no_personal_or_machine_identifying_data(tmp_path):
    import getpass
    import os
    import socket

    path = run(load_experiment(EXAMPLE)).manifest.write(tmp_path / "manifest.json")
    text = path.read_text()

    leaks = {
        "hostname": socket.gethostname(),
        "username": getpass.getuser(),
        "home directory": str(Path.home()),
        "working directory": os.getcwd(),
    }
    for description, value in leaks.items():
        if value and len(value) > 3:
            assert value not in text, f"manifest leaked the {description}"

    manifest = json.loads(text)
    assert set(manifest["environment"]) == {
        "python",
        "implementation",
        "os_family",
        "numpy",
        "marse",
    }
    assert manifest["started_utc"].endswith("+00:00")
    assert manifest["finished_utc"].endswith("+00:00")
    # No absolute paths anywhere in the manifest.
    assert "/home/" not in text
    assert "C:\\" not in text


def test_manifest_records_the_models_and_streams_used(tmp_path):
    manifest = run(load_experiment(EXAMPLE)).manifest
    assert manifest.models["growth"] == "monod_v1"
    assert manifest.models["integrator"] == "rk4_v1"
    assert "core" in manifest.random_streams
    assert manifest.seed == 20260924


def test_manifest_json_is_stable_across_writes(tmp_path):
    """Byte-identical output for the same run, so manifests diff cleanly."""
    config = load_experiment(EXAMPLE)
    manifest = run(config).manifest
    first = manifest.write(tmp_path / "a.json").read_text()
    second = manifest.write(tmp_path / "b.json").read_text()
    assert first == second


# --- the command-line interface ----------------------------------------------


def marse(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "marse", *args],
        capture_output=True,
        text=True,
        timeout=600,
        cwd=cwd,
    )


def test_cli_run_then_replay(tmp_path):
    outputs = tmp_path / "out"
    run_result = marse("run", str(EXAMPLE), "-o", str(outputs))
    assert run_result.returncode == 0, run_result.stderr
    assert (outputs / "manifest.json").is_file()
    assert (outputs / "trajectory.csv").is_file()

    replay = marse("replay", str(outputs / "manifest.json"))
    assert replay.returncode == 0, replay.stderr
    assert "reproduced the recorded final state exactly" in replay.stdout


def test_cli_trajectory_has_units_in_its_headers(tmp_path):
    outputs = tmp_path / "out"
    marse("run", str(EXAMPLE), "-o", str(outputs))
    header = (outputs / "trajectory.csv").read_text().splitlines()[0]
    assert header == "time_h,glucose_mm,fast_grower_g_per_l,gleaner_g_per_l"


def test_cli_reports_an_invalid_experiment_without_a_traceback(tmp_path):
    broken = tmp_path / "broken.json"
    broken.write_text(json.dumps(minimal_experiment(duration_h=-1.0)))
    result = marse("run", str(broken))
    assert result.returncode == 2
    assert "invalid experiment" in result.stdout
    assert "Traceback" not in result.stderr


def test_cli_replay_detects_a_tampered_manifest(tmp_path):
    outputs = tmp_path / "out"
    marse("run", str(EXAMPLE), "-o", str(outputs))
    path = outputs / "manifest.json"
    raw = json.loads(path.read_text())
    raw["config"]["organisms"][0]["mu_opt_per_h"] = 2.0
    path.write_text(json.dumps(raw))

    result = marse("replay", str(path))
    assert result.returncode == 2
    assert "checksum" in result.stdout


def test_cli_with_no_command_prints_help():
    result = marse()
    assert result.returncode == 0
    assert "run a simulation" in result.stdout


def test_config_is_a_frozen_dataclass():
    config = load_experiment(EXAMPLE)
    assert isinstance(config, ExperimentConfig)
    with pytest.raises((AttributeError, TypeError)):
        config.duration_h = 1.0  # type: ignore[misc]
