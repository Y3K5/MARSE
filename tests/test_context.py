"""Source-scoped context contracts; fixtures are entirely synthetic."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import pytest

from marse.context import ContextError, validate_context_exchange

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def bundle():
    return json.loads((ROOT / "examples/context_exchange/synthetic.json").read_text())


def test_synthetic_context_is_valid_and_input_is_unchanged(bundle):
    before = copy.deepcopy(bundle)
    summary = validate_context_exchange(bundle)
    assert summary.contexts == 1
    assert summary.organisms == 2
    assert summary.parameters == 1
    assert summary.relationships == 1
    assert summary.molecular_execution == "not_implemented"
    assert bundle == before


@pytest.mark.parametrize(
    ("collection", "field", "value", "error"),
    [
        ("parameters", "units", None, "Missing units"),
        ("parameters", "value", 0, "Unknown must remain null"),
        ("parameters", "context_id", "missing", "Unknown contexts"),
        ("parameters", "host", "canine", "Host/site/state"),
        ("parameters", "site", "other", "Host/site/state"),
        ("parameters", "disease_context", "healthy", "Host/site/state"),
        ("parameters", "compartment", "other", "Unknown compartment"),
        ("parameters", "evidence_class", "measured", "parameter value"),
        ("relationships", "native_access_claim", True, "Native access"),
        ("relationships", "field_overlap_is_contact", True, "Overlap is not contact"),
        ("relationships", "coupled_to_solver", True, "coupling not implemented"),
        ("relationships", "structure_status", "resolved", "Structure link"),
        ("relationships", "mediator_id", "invented", "Mediator mapping"),
        ("relationships", "participants", ["anonymous-a", "missing"], "Unknown organisms"),
        ("relationships", "participants", ["anonymous-a", "anonymous-a"], "Duplicate participant"),
        ("contexts", "biofilm_age", None, "biofilm_age"),
        ("contexts", "disease_context", "mature", "Unknown disease"),
        ("evidence", "host", "canine", "Host transfer"),
        ("evidence", "observation", "modified", "Annotation checksum"),
        ("organisms", "host", "canine", "Host transfer"),
    ],
)
def test_incompatible_records_are_rejected(bundle, collection, field, value, error):
    bundle[collection][0][field] = value
    with pytest.raises(ContextError, match=error):
        validate_context_exchange(bundle)


def test_duplicate_ids_rejected(bundle):
    bundle["organisms"].append(copy.deepcopy(bundle["organisms"][0]))
    with pytest.raises(ContextError, match="Duplicate organisms ID"):
        validate_context_exchange(bundle)


def test_molecular_transform_cannot_be_invented(bundle):
    bundle["contexts"][0]["frame"]["molecular_transform"] = [1, 0, 0, 1]
    with pytest.raises(ContextError, match="Molecular transform"):
        validate_context_exchange(bundle)


def test_fasta_file_digest_cannot_be_a_normalized_sequence_digest(bundle):
    organism = bundle["organisms"][0]
    organism["accession"] = "SYNTHETIC"
    organism["sequence_digest"] = {"kind": "fasta_file_sha256", "value": "a" * 64}
    with pytest.raises(ContextError, match="Digest type"):
        validate_context_exchange(bundle)


def test_annotation_digest_cannot_be_a_source_file_digest(bundle):
    bundle["evidence"][0]["source_digest"] = bundle["evidence"][0]["annotation_digest"]
    with pytest.raises(ContextError, match="Digest type"):
        validate_context_exchange(bundle)


def test_execution_cannot_claim_outputs_without_a_digest(bundle):
    bundle["contexts"][0]["statuses"]["execution"] = "executed_synthetic"
    with pytest.raises(ContextError, match="Missing executed output"):
        validate_context_exchange(bundle)


def test_all_five_status_axes_are_required(bundle):
    del bundle["contexts"][0]["statuses"]["evidence"]
    with pytest.raises(ContextError, match="Five independent"):
        validate_context_exchange(bundle)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), True])
def test_nonfinite_or_boolean_measurements_rejected(bundle, value):
    parameter = bundle["parameters"][0]
    parameter["evidence_class"] = "hypothetical_scenario"
    parameter["value"] = value
    with pytest.raises(ContextError, match="parameter value"):
        validate_context_exchange(bundle)


def test_observation_compartment_cannot_be_silently_transferred(bundle):
    parameter = bundle["parameters"][0]
    parameter.update(evidence_class="measured", value=1.0, evidence_ids=["synthetic-observation"])
    bundle["evidence"][0]["site"] = "saliva"
    with pytest.raises(ContextError, match="Evidence compartment"):
        validate_context_exchange(bundle)


def test_context_can_keep_explicit_scenario_parameters(bundle):
    bundle["parameters"][0].update(evidence_class="hypothetical_scenario", value=1.0)
    assert validate_context_exchange(bundle).parameters == 1


def test_synthetic_example_hashes_bind_actual_reference_and_protocol(bundle):
    engine = bundle["contexts"][0]["engine"]
    for key, name in (("source_digest", "reference.cjs"), ("config_digest", "protocol.json")):
        actual = hashlib.sha256((ROOT / "examples/numerical_reference" / name).read_bytes())
        assert engine[key]["value"] == actual.hexdigest()


def test_unknown_schema_is_rejected(bundle):
    bundle["schema"] = "future/2"
    with pytest.raises(ContextError, match="Unsupported context"):
        validate_context_exchange(bundle)


@pytest.mark.parametrize(
    ("collection", "field", "value"),
    [
        ("evidence", "host", []),
        ("organisms", "host", {}),
        ("contexts", "host", []),
        ("contexts", "disease_context", []),
        ("parameters", "evidence_class", []),
        ("relationships", "evidence_class", {}),
    ],
)
def test_malformed_choice_types_raise_contract_errors(bundle, collection, field, value):
    bundle[collection][0][field] = value
    with pytest.raises(ContextError):
        validate_context_exchange(bundle)
