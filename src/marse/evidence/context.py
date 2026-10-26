"""Validate evidence context without executing a scientific model.

This narrow exchange keeps source scope, units and missing measurements explicit.
It does not load molecular coordinates, infer contact or run a host-response model.
Only standard-library code is used; inputs are inspected without mutation.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

__all__ = ["ContextError", "ContextSummary", "validate_context_exchange"]

SCHEMA = "marse.context-exchange/1"
COLLECTIONS = ("evidence", "contexts", "organisms", "parameters", "relationships")
HOSTS = {"synthetic", "human", "canine"}
STATUS_VALUES = {
    "implementation": {"recorded", "implemented", "not_implemented"},
    "execution": {"not_run", "recorded", "executed_synthetic", "saved_legacy"},
    "gate": {"not_evaluated", "passed", "failed", "unresolved"},
    "evidence": {"synthetic_only", "scenario_only", "source_scoped", "unresolved"},
    "approval": {"local_only", "review_required"},
}


class ContextError(ValueError):
    """An exchange has malformed data or an unsupported inference."""


@dataclass(frozen=True, slots=True)
class ContextSummary:
    """Counts of validated records, with no biological qualification implied."""

    contexts: int
    evidence: int
    organisms: int
    parameters: int
    relationships: int
    molecular_execution: str = "not_implemented"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContextError(message)


def _object(value: Any, label: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), f"{label} must be an object")
    return value


def _text(value: Any, label: str) -> str:
    _require(isinstance(value, str) and bool(value.strip()), f"Missing {label}")
    return value


def _number(value: Any, label: str, minimum: float = 0.0) -> None:
    _require(
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(value)
        and value >= minimum,
        f"Invalid {label}",
    )


def _digest(value: Any, kind: str) -> None:
    record = _object(value, "digest")
    _require(record.get("kind") == kind, "Digest type mismatch")
    _require(
        isinstance(record.get("value"), str)
        and re.fullmatch(r"[0-9a-f]{64}", record["value"]) is not None,
        "Invalid SHA-256",
    )


def _statuses(record: Mapping[str, Any]) -> None:
    statuses = _object(record.get("statuses"), "statuses")
    _require(set(statuses) == set(STATUS_VALUES), "Five independent statuses required")
    for key, choices in STATUS_VALUES.items():
        _require(isinstance(statuses[key], str) and statuses[key] in choices, f"Invalid {key}")


def validate_context_exchange(bundle: Mapping[str, Any]) -> ContextSummary:
    """Check identities, provenance, units, host/site scope and unresolved links.

    Parameters
    ----------
    bundle
        A decoded ``marse.context-exchange/1`` document. See the synthetic example
        and ``docs/context-exchange.md`` for the required record fields.

    Returns
    -------
    ContextSummary
        Record counts. Validation establishes contract consistency only.

    Raises
    ------
    ContextError
        If required fields, joins or claim boundaries are invalid.
    """
    bundle = _object(bundle, "exchange")
    _require(bundle.get("schema") == SCHEMA, "Unsupported context schema")
    _require(set(bundle) == {"schema", *COLLECTIONS}, "Unexpected exchange collections")
    maps: dict[str, dict[str, Mapping[str, Any]]] = {}
    for collection in COLLECTIONS:
        records = bundle.get(collection)
        _require(isinstance(records, list), f"Missing {collection}")
        maps[collection] = {}
        for record in records:
            record = _object(record, collection)
            identity = _text(record.get("id"), "record ID")
            _require(re.fullmatch(r"[a-z0-9][a-z0-9_.-]*", identity) is not None, "Invalid ID")
            _require(identity not in maps[collection], f"Duplicate {collection} ID")
            _statuses(record)
            maps[collection][identity] = record

    def reference(collection: str, identity: Any) -> Mapping[str, Any]:
        _text(identity, "reference ID")
        _require(identity in maps[collection], f"Unknown {collection} ID")
        return maps[collection][identity]

    def scope(record: Mapping[str, Any], context: Mapping[str, Any]) -> None:
        _require(
            all(record.get(key) == context[key] for key in ("host", "site", "disease_context")),
            "Host/site/state mismatch",
        )

    for evidence in maps["evidence"].values():
        _require(
            isinstance(evidence.get("host"), str) and evidence["host"] in HOSTS,
            "Unknown evidence host",
        )
        for key in ("source", "site", "method", "observation", "limits", "independence_group"):
            _text(evidence.get(key), key)
        annotation = evidence.get("annotation_digest")
        _digest(annotation, "annotation_text_sha256")
        _require(
            annotation["value"] == hashlib.sha256(evidence["observation"].encode()).hexdigest(),
            "Annotation checksum mismatch",
        )
        source_digest = evidence.get("source_digest")
        if source_digest is not None:
            _digest(source_digest, "source_file_sha256")

    for organism in maps["organisms"].values():
        _require(
            isinstance(organism.get("host"), str) and organism["host"] in HOSTS,
            "Unknown organism host",
        )
        for key in ("name", "strain", "representation"):
            _text(organism.get(key), key)
        sequence_digest = organism.get("sequence_digest")
        if sequence_digest is not None:
            _text(organism.get("accession"), "accession")
            _digest(sequence_digest, "normalized_sequence_sha256")

    for context in maps["contexts"].values():
        _require(
            isinstance(context.get("host"), str) and context["host"] in HOSTS,
            "Unknown context host",
        )
        _require(
            isinstance(context.get("disease_context"), str)
            and context["disease_context"]
            in {"healthy", "gingivitis", "periodontitis", "unresolved"},
            "Unknown disease context",
        )
        for key in ("site", "biofilm_age", "biofilm_phase", "physiology", "host_response"):
            _text(context.get(key), key)
        _require(context.get("representation") == "density_fields_not_cells", "Invalid scale")
        frame = _object(context.get("frame"), "frame")
        _text(frame.get("id"), "frame ID")
        _text(frame.get("units"), "frame units")
        _require(
            "molecular_transform" in frame and frame["molecular_transform"] is None,
            "Molecular transform not implemented; explicit null required",
        )
        engine = _object(context.get("engine"), "engine")
        for key in ("id", "version", "run_id", "time_units", "seed_policy"):
            _text(engine.get(key), key)
        for key, kind in (
            ("source_digest", "source_file_sha256"),
            ("config_digest", "config_file_sha256"),
        ):
            _digest(engine.get(key), kind)
        if engine.get("output_digest") is None:
            _require(
                context["statuses"]["execution"] == "not_run", "Missing executed output digest"
            )
        else:
            _digest(engine["output_digest"], "output_file_sha256")
        organisms = context.get("organism_ids")
        _require(isinstance(organisms, list), "Missing organism IDs")
        for identity in organisms:
            _text(identity, "organism ID")
        _require(len(organisms) == len(set(organisms)), "Duplicate context organism")
        for identity in organisms:
            _require(reference("organisms", identity)["host"] == context["host"], "Host transfer")

    for parameter in maps["parameters"].values():
        context = reference("contexts", parameter.get("context_id"))
        scope(parameter, context)
        for key in (
            "compartment",
            "units",
            "method",
            "uncertainty",
            "calibration_role",
            "identifiability",
            "transfer_limits",
            "needed_observation",
        ):
            _text(parameter.get(key), key)
        _require(
            parameter["compartment"]
            in {
                "saliva",
                "entrance",
                "lumen",
                "plaque_interior",
                "epithelium",
                "tissue",
            },
            "Unknown compartment",
        )
        evidence_class = parameter.get("evidence_class")
        _require(
            isinstance(evidence_class, str)
            and evidence_class in {"measured", "inferred", "hypothetical_scenario", "unresolved"},
            "Unknown parameter evidence class",
        )
        _require("value" in parameter, "Missing parameter value or explicit null")
        if evidence_class == "unresolved":
            _require(parameter["value"] is None, "Unknown must remain null")
        else:
            _number(parameter["value"], "parameter value", minimum=-math.inf)
        evidence_ids = parameter.get("evidence_ids")
        _require(isinstance(evidence_ids, list), "Missing parameter evidence IDs")
        if evidence_class in {"measured", "inferred"}:
            _require(bool(evidence_ids), "Measured/inferred value needs source evidence")
        for identity in evidence_ids:
            evidence = reference("evidence", identity)
            _require(evidence["host"] == context["host"], "Evidence host transfer")
            _require(evidence["site"] == parameter["compartment"], "Evidence compartment mismatch")

    for relationship in maps["relationships"].values():
        context = reference("contexts", relationship.get("context_id"))
        scope(relationship, context)
        participants = relationship.get("participants")
        _require(isinstance(participants, list) and len(participants) >= 2, "Missing participants")
        for identity in participants:
            reference("organisms", identity)
            _require(identity in context["organism_ids"], "Participant outside context")
        _require(len(set(participants)) == len(participants), "Duplicate participant")
        _require(
            isinstance(relationship.get("evidence_class"), str)
            and relationship["evidence_class"]
            in {
                "co_occurrence",
                "spatial_neighbors",
                "functional_perturbation",
                "resource_exchange",
                "host_mediated",
                "unresolved",
            },
            "Unknown relationship evidence class",
        )
        for key in ("conditions", "direction", "alternatives", "limits"):
            _text(relationship.get(key), key)
        evidence_ids = relationship.get("evidence_ids")
        _require(isinstance(evidence_ids, list), "Missing relationship evidence IDs")
        for identity in evidence_ids:
            _require(reference("evidence", identity)["host"] == context["host"], "Host transfer")
        _require(
            relationship.get("coupled_to_solver") is False, "Relationship coupling not implemented"
        )
        _require(relationship.get("field_overlap_is_contact") is False, "Overlap is not contact")
        _require(
            relationship.get("native_access_claim") is False, "Native access is not established"
        )
        _require(relationship.get("mediator_id") is None, "Mediator mapping not implemented")
        _require(relationship.get("structure_status") == "unresolved", "Structure link unresolved")

    return ContextSummary(**{key: len(value) for key, value in maps.items()})
