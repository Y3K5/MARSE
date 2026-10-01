"""Ecosystem runs are recorded and replayed like every other MARSE run.

A manifest must rebuild exactly the configuration that ran, so the first tests
here prove that an ecosystem configuration survives the trip through JSON
with every feature switched on, defaults included. The rest check what the
manifest records and that a replay reproduces the run bit for bit.
"""

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from marse.cli import main
from marse.core.config import load_experiment as load_core_experiment
from marse.core.provenance import Manifest, config_checksum
from marse.core.simulation import run as run_core
from marse.ecosystem import (
    AdditiveConfig,
    ConditionConfig,
    EcosystemConfig,
    NutrientConfig,
    PhenotypeConfig,
    SeedRegion,
    SpeciesConfig,
    run,
)
from marse.ecosystem.model import ENGINE_VERSION, _state_digest, ecosystem_from_dict
from marse.experimental.host.immune import ImmuneInteraction, MolecularNeutralizer
from marse.microbes.additives import AdditiveEffect
from marse.microbes.niche import Capability

EXPERIMENTS = Path(__file__).resolve().parents[1] / "examples" / "experiments"


def everything_on() -> EcosystemConfig:
    """A small configuration with every optional feature set to a non-default."""
    capability = Capability(
        "carbon-use",
        0.6,
        "carbon",
        0.3,
        temperature_c=(5.0, 37.0, 45.0),
        ph=(4.0, 7.0, 9.0),
        oxygen_half_saturation=0.05,
        evidence_source="test",
        evidence_confidence="C",
    )
    return EcosystemConfig(
        experiment_id="everything-on",
        width=5,
        height=4,
        cell_size_um=10.0,
        duration_h=0.3,
        timestep_h=0.01,
        seed=11,
        nutrients=(
            NutrientConfig("oxygen", 0.2, 10.0, 0.25, ("top",)),
            NutrientConfig("carbon", 1.0, 5.0),
        ),
        conditions=(
            ConditionConfig("temperature_c", 35.0, 0.0),
            ConditionConfig("ph", 7.1, 0.0),
        ),
        additives=(
            AdditiveConfig("signal", 0.0, 2.0, decay_per_h=0.1),
            AdditiveConfig("effector", 0.5, 1.0, 0.0, 1.0, ("bottom",)),
            AdditiveConfig("decoy", 0.2, 1.0),
        ),
        species=(
            SpeciesConfig(
                "builder",
                0.05,
                0.7,
                (0.1, 0.3),
                (1.0, 2.0),
                mutation_probability=0.2,
                mutation_growth_multiplier=1.1,
                seed_regions=(SeedRegion(0.2, 0.8, 0.3, 0.4),),
                spreading_per_h=5.0,
                production_per_nutrient=(0.0, 0.01),
                capabilities=(capability,),
                additive_effects=(AdditiveEffect("effector", 0.2, 1.0, 0.5, 2.0, "decreasing"),),
                chemotaxis_field="signal",
                chemotaxis_sensitivity=0.5,
                phenotypes=(PhenotypeConfig("dense", 0.3, 0.1, 0.8, 0.5, 0.05),),
                quorum_signal="signal",
                quorum_secretion_per_h=0.2,
                adhesion_per_h=0.1,
                detachment_per_h=0.05,
                adhesion_edges=("bottom",),
            ),
            SpeciesConfig("plain", 0.02, 0.4, (0.2, 0.2), (1.0, 1.0)),
        ),
        mutation_interval_h=0.1,
        carrying_capacity=2.0,
        competition_coefficients=((0.5, 0.2), (0.3, 0.5)),
        immune_interactions=(ImmuneInteraction("plain", "effector", 0.5, 0.4, 1.5, 0.8),),
        immune_neutralizers=(MolecularNeutralizer("effector", "decoy", 0.5, 0.3, 2.0),),
    )


def through_json(config: EcosystemConfig) -> EcosystemConfig:
    return ecosystem_from_dict(json.loads(json.dumps(config.to_dict())))


@pytest.mark.parametrize(
    "name", ["two_species_ecosystem.json", "periodontal_pathogen_biofilm.json"]
)
def test_example_configurations_survive_the_trip_through_json(name: str):
    raw = json.loads((EXPERIMENTS / name).read_text("utf-8"))
    config = ecosystem_from_dict(raw)
    again = through_json(config)
    assert config_checksum(again) == config_checksum(config)
    assert again.to_dict() == config.to_dict()


def test_every_optional_feature_survives_the_trip_through_json():
    config = everything_on()
    again = through_json(config)
    # Compared as records, not objects: an omitted production list comes back
    # as the explicit zeros it meant, which is the point of writing it out.
    assert again.to_dict() == config.to_dict()
    assert config_checksum(again) == config_checksum(config)
    assert again.species[0] == config.species[0]  # every optional field set, so exact


def test_implied_defaults_are_written_out():
    # An omitted production list means zeros; the manifest records the zeros.
    species = everything_on().to_dict()["species"][1]
    assert species["production_per_nutrient"] == [0.0, 0.0]
    assert species["mutation_probability"] == 0.0


def test_the_manifest_records_what_ran():
    result = run(everything_on())
    manifest = result.manifest
    assert manifest.kind == "ecosystem"
    assert manifest.models["engine"] == ENGINE_VERSION
    assert set(result.provider_versions.items()) <= set(manifest.models.items())
    assert manifest.random_streams == ("ecosystem.mutations",)
    assert manifest.steps == everything_on().steps
    assert set(manifest.outputs) == {
        "final_time_h",
        "biomass_summed_over_cells",
        "nutrient_summed_over_cells",
        "mutated_cells",
        "final_state_sha256",
    }
    assert manifest.outputs["mutated_cells"]["builder"] > 0


def test_a_replay_through_the_manifest_is_bit_identical(tmp_path: Path):
    original = run(everything_on())
    path = original.manifest.write(tmp_path / "manifest.json")
    replayed = run(Manifest.read(path).experiment())
    assert replayed.manifest.run_id == original.manifest.run_id
    assert replayed.manifest.outputs == original.manifest.outputs
    for field in ("biomass", "nutrients", "additives", "mutations", "phenotype_indices"):
        np.testing.assert_array_equal(
            getattr(replayed.final_state, field), getattr(original.final_state, field)
        )


def test_the_state_digest_sees_what_the_totals_miss():
    # Mirror the biomass left to right: every total is unchanged, the state is not.
    state = run(everything_on()).final_state
    mirrored = replace(state, biomass=state.biomass[:, :, ::-1])
    assert np.allclose(mirrored.biomass.sum(axis=(1, 2)), state.biomass.sum(axis=(1, 2)))
    assert _state_digest(mirrored) != _state_digest(state)
    assert _state_digest(replace(state)) == _state_digest(state)


def test_a_different_seed_is_a_different_run():
    first = run(everything_on()).manifest
    second = run(replace(everything_on(), seed=12)).manifest
    assert second.run_id != first.run_id
    assert second.outputs["final_state_sha256"] != first.outputs["final_state_sha256"]


def test_an_edited_manifest_is_refused(tmp_path: Path):
    path = run(everything_on()).manifest.write(tmp_path / "manifest.json")
    raw = json.loads(path.read_text("utf-8"))
    raw["config"]["carrying_capacity"] = 3.0
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        Manifest.read(path).experiment()


def test_a_manifest_whose_kind_contradicts_its_configuration_is_refused(tmp_path: Path):
    # A batch configuration labelled as a biofilm run parses, so only the kind
    # check can catch it.
    result = run_core(load_core_experiment(EXPERIMENTS / "two_species_batch.json"))
    raw = result.manifest.to_dict()
    raw["kind"] = "biofilm_profile"
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="declares a 'biofilm_profile' run"):
        Manifest.read(path).experiment()


def test_an_ecosystem_configuration_is_not_replayed_by_the_batch_engine(tmp_path: Path):
    path = run(everything_on()).manifest.write(tmp_path / "manifest.json")
    raw = json.loads(path.read_text("utf-8"))
    raw["kind"] = "batch"
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown field"):
        Manifest.read(path).experiment()


def test_an_unknown_kind_is_refused(tmp_path: Path):
    path = run(everything_on()).manifest.write(tmp_path / "manifest.json")
    raw = json.loads(path.read_text("utf-8"))
    raw["kind"] = "agent_based"
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="kind 'agent_based' is not supported"):
        Manifest.read(path).experiment()


def test_version_one_manifests_are_still_readable(tmp_path: Path):
    for experiment, kind in (
        ("two_species_batch.json", "batch"),
        ("biofilm_oxygen_profile.json", "biofilm_profile"),
    ):
        result = run_core(load_core_experiment(EXPERIMENTS / experiment))
        raw = result.manifest.to_dict()
        raw["manifest_version"] = 1
        del raw["kind"]
        path = tmp_path / f"{kind}.json"
        path.write_text(json.dumps(raw), encoding="utf-8")
        old = Manifest.read(path)
        assert old.kind == kind
        assert old.experiment().kind == kind


def test_the_manifest_holds_no_path_or_machine_identity(tmp_path: Path):
    text = run(everything_on()).manifest.write(tmp_path / "manifest.json").read_text("utf-8")
    assert str(tmp_path) not in text
    assert str(Path.home()) not in text


def test_cli_replay_reports_a_recorded_result_that_no_longer_matches(tmp_path: Path, capsys):
    path = tmp_path / "experiment.json"
    path.write_text(json.dumps(everything_on().to_dict()), encoding="utf-8")
    assert main(["ecosystem", str(path), "-o", str(tmp_path / "run")]) == 0
    manifest = tmp_path / "run" / "manifest.json"
    assert main(["replay", str(manifest)]) == 0
    assert "reproduced the recorded results exactly" in capsys.readouterr().out

    # Outputs are not covered by the configuration checksum, so an edited
    # result reads cleanly - and the replay has to catch it.
    raw = json.loads(manifest.read_text("utf-8"))
    raw["outputs"]["biomass_summed_over_cells"]["plain"] += 1.0
    manifest.write_text(json.dumps(raw), encoding="utf-8")
    assert main(["replay", str(manifest)]) == 1
    assert "DIFFERS" in capsys.readouterr().out
