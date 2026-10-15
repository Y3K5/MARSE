"""The JSON feed an importer writes, and how it meets an experiment file."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from marse.ingest.errors import IngestError

FORMAT = "marse.ingest/1"


@dataclass(frozen=True, slots=True)
class Feed:
    """A checked table. ``experiment_fragment`` is null when nothing maps."""

    kind: str
    source_name: str
    columns: tuple[str, ...]
    rows: tuple[dict[str, str], ...]
    experiment_fragment: dict[str, Any] | None
    limits: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "columns": list(self.columns),
            "experiment_fragment": self.experiment_fragment,
            "format": FORMAT,
            "kind": self.kind,
            "limits": self.limits,
            "qc": {"checks": _checks(self.kind), "passed": True, "rows": len(self.rows)},
            "records": [dict(row) for row in self.rows],
            "source_name": self.source_name,
        }


def _checks(kind: str) -> list[str]:
    shared = ["header present", "values finite"]
    if kind == "environment":
        return [*shared, "temperature_c in -20..150 C", "ph in 0..14", "units C and pH"]
    if kind == "abundance":
        return [*shared, "abundance non-negative", "relative rows sum to 1 per sample"]
    if kind == "transcriptomics":
        return [*shared, "value non-negative", "unit count or tpm"]
    return [*shared, "value non-negative", "unit mmol_per_l, mol_per_m3, or arbitrary"]


def apply_fragment(experiment: dict[str, Any], fragment: dict[str, Any]) -> dict[str, Any]:
    """Copy fragment keys onto a new experiment dict.

    A key that is already present and differs is refused. The call does not
    invent organisms, yields, or a domain.
    """
    if not isinstance(experiment, dict) or not isinstance(fragment, dict):
        raise IngestError("an experiment and a fragment must both be mappings")
    merged = dict(experiment)
    for key, value in fragment.items():
        if key in merged and merged[key] != value:
            raise IngestError(f"experiment already has {key}")
        merged[key] = value
    return merged
