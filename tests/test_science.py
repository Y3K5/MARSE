"""Tests for evidence-backed culture records."""

from dataclasses import replace

import pytest

from marse.science import (
    AgarProtocol,
    CultureDataset,
    CultureError,
    CultureRecord,
    KineticObservation,
    MeasurementMethod,
    MediumRecipe,
    SourceRecord,
    compile_culture,
    dataset_checksum,
    load_culture_dataset,
)


def dataset() -> CultureDataset:
    source = SourceRecord(
        id="method-standard",
        title="Aerobic plate count",
        organization_or_authors="FDA",
        year=2024,
        source_type="method_standard",
        evidence_grade="B",
        url="https://www.fda.gov/food/laboratory-methods-food/bam-chapter-3-aerobic-plate-count",
    )
    medium = MediumRecipe(
        id="defined-medium",
        name="Defined carbon medium",
        components={"glucose": "4 g/L", "nitrogen": "1 g/L"},
        ph=7.0,
        source_ids=("method-standard",),
    )
    measurement = MeasurementMethod(
        id="colony-count",
        endpoint="cfu_ml",
        unit="CFU/mL",
        instrument_or_method="serial dilution plate count",
    )
    culture = CultureRecord(
        id="strain-solid-assay",
        organism="example bacterium",
        strain="strain-1",
        medium_id="defined-medium",
        format="agar",
        temperature_c=37.0,
        ph=7.0,
        atmosphere="aerobic",
        inoculum="overnight culture, 1:100 dilution",
        incubation_h=24.0,
        measurement_id="colony-count",
        source_ids=("method-standard",),
        agar=AgarProtocol(15.0, 4.0, 20.0, "surface spread", "aerobic"),
        kinetics=KineticObservation(
            measurement_id="colony-count",
            model="logistic",
            growth_rate_per_h=0.8,
            lag_h=2.0,
            carrying_capacity=1.0e9,
            uncertainty=0.1,
            replicates=3,
            fitting_interval_h=(0.0, 24.0),
        ),
    )
    return CultureDataset("2026.1", (source,), (medium,), (measurement,), (culture,))


def test_dataset_validates_references_and_has_stable_checksum(tmp_path):
    evidence = dataset()
    assert dataset_checksum(evidence) == dataset_checksum(evidence)
    path = evidence.write_json(tmp_path / "culture.json")
    loaded = load_culture_dataset(path)
    assert loaded == evidence
    assert dataset_checksum(loaded) == dataset_checksum(evidence)


def test_agar_requires_protocol_and_non_universal_context():
    with pytest.raises(CultureError, match="agar protocol"):
        CultureRecord(
            id="missing-agar",
            organism="bacterium",
            strain="strain-1",
            medium_id="medium",
            format="agar",
            temperature_c=37.0,
            ph=7.0,
            atmosphere="aerobic",
            inoculum="fresh culture",
            incubation_h=24.0,
            measurement_id="count",
            source_ids=("source",),
        )

    with pytest.raises(CultureError, match="universal organism"):
        CultureRecord(
            id="universal",
            organism="all",
            strain=None,
            medium_id="medium",
            format="broth",
            temperature_c=37.0,
            ph=7.0,
            atmosphere="aerobic",
            inoculum="fresh culture",
            incubation_h=24.0,
            measurement_id="count",
            source_ids=("source",),
        )


def test_optimum_claim_requires_strain_and_od_requires_calibration():
    with pytest.raises(CultureError, match="optimum claim"):
        CultureRecord(
            id="unscoped-optimum",
            organism="bacterium",
            strain=None,
            medium_id="medium",
            format="broth",
            temperature_c=37.0,
            ph=7.0,
            atmosphere="aerobic",
            inoculum="fresh culture",
            incubation_h=24.0,
            measurement_id="count",
            source_ids=("source",),
            optimum_claim=True,
        )

    with pytest.raises(CultureError, match="calibration"):
        MeasurementMethod("od", "od", "AU", "plate reader")


def test_dataset_rejects_unknown_source_and_medium():
    with pytest.raises(CultureError, match="unknown source"):
        CultureDataset(
            "2026.1",
            (),
            (
                MediumRecipe(
                    "medium",
                    "Medium",
                    {"carbon": "1 g/L"},
                    source_ids=("missing",),
                ),
            ),
            (),
            (),
        )


def test_compiler_creates_a_traceable_batch_experiment():
    evidence = dataset()
    culture = replace(evidence.cultures[0], format="broth", agar=None)
    evidence = replace(evidence, cultures=(culture,))
    compilation = compile_culture(
        evidence,
        culture.id,
        substrate_name="glucose",
        substrate_initial_mm=10.0,
        initial_biomass_g_per_l=0.01,
        yield_g_per_mmol=0.1,
        k_s_mm=0.2,
        duration_h=8.0,
        timestep_h=0.1,
        seed=7,
    )
    assert compilation.runnable
    assert compilation.experiment is not None
    assert compilation.experiment.organisms[0].mu_opt_per_h == 0.8
    assert any("evidence dataset sha256=" in item for item in compilation.assumptions)
    assert compilation.selected_source_ids == ("method-standard",)


def test_compiler_refuses_to_project_agar_into_well_mixed_batch():
    compilation = compile_culture(
        dataset(),
        "strain-solid-assay",
        substrate_name="glucose",
        substrate_initial_mm=10.0,
        initial_biomass_g_per_l=0.01,
        yield_g_per_mmol=0.1,
        k_s_mm=0.2,
        duration_h=8.0,
        timestep_h=0.1,
    )
    assert not compilation.runnable
    assert compilation.experiment is None
    assert "spatial compiler" in compilation.unresolved[0]


def test_documented_culture_evidence_dataset_loads():
    path = "examples/experiments/culture_evidence.json"
    evidence = load_culture_dataset(path)
    assert evidence.dataset_version == "2026.1"
    assert len(evidence.sources) == 2
    assert evidence.cultures == ()
