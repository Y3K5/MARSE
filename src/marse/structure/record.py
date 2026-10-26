"""A predicted structure and the scores that describe it."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any


class StructureError(ValueError):
    """A sequence or a score payload cannot be checked."""


_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True, slots=True)
class StructureQC:
    """Range checks and the summary AlphaFold already publishes.

    ``confident_fraction`` uses pLDDT at least 70, and ``high_fraction`` uses
    at least 90. Those cuts are AlphaFold's reported bins, not a MARSE result.
    """

    mean_plddt: float
    min_plddt: float
    high_fraction: float
    confident_fraction: float
    pae_present: bool
    mean_pae: float | None
    max_pae: float | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "confident_fraction": self.confident_fraction,
            "high_fraction": self.high_fraction,
            "max_pae": self.max_pae,
            "mean_pae": self.mean_pae,
            "mean_plddt": self.mean_plddt,
            "min_plddt": self.min_plddt,
            "pae_present": self.pae_present,
        }


@dataclass(frozen=True, slots=True)
class StructureRecord:
    """One prediction. ``limits`` says what the scores do not decide."""

    adapter: str
    sequence: str
    plddt: tuple[float, ...]
    pae: tuple[tuple[float, ...], ...] | None
    qc: StructureQC
    limits: str = "pLDDT and PAE do not set a growth rate, a yield, or a biomass formula."

    def to_dict(self) -> dict[str, Any]:
        return {
            "adapter": self.adapter,
            "length": len(self.sequence),
            "limits": self.limits,
            "pae_present": self.qc.pae_present,
            "qc": self.qc.to_dict(),
        }


def check_sequence(sequence: str) -> str:
    """Accept the twenty standard residues. Ambiguity codes are refused."""
    if not sequence or any(residue not in _AMINO_ACIDS for residue in sequence):
        raise StructureError("sequence must be standard amino acids, with no ambiguity codes")
    return sequence


def score_prediction(
    adapter: str,
    sequence: str,
    plddt: list[Any] | tuple[Any, ...],
    pae: list[Any] | None,
) -> StructureRecord:
    """Build a record when every supplied score is finite and in range."""
    sequence = check_sequence(sequence)
    if len(plddt) != len(sequence):
        raise StructureError("plddt must have one value per residue")
    scores = tuple(_plddt(value, index) for index, value in enumerate(plddt))
    matrix = None if pae is None else _pae(pae, len(sequence))
    count = len(scores)
    return StructureRecord(
        adapter,
        sequence,
        scores,
        matrix,
        StructureQC(
            mean_plddt=sum(scores) / count,
            min_plddt=min(scores),
            high_fraction=sum(value >= 90 for value in scores) / count,
            confident_fraction=sum(value >= 70 for value in scores) / count,
            pae_present=matrix is not None,
            mean_pae=None if matrix is None else _mean_pae(matrix),
            max_pae=None if matrix is None else max(max(row) for row in matrix),
        ),
    )


def _plddt(value: float, index: int) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise StructureError(f"plddt[{index}] must be a number")
    number = float(value)
    if not math.isfinite(number) or not 0.0 <= number <= 100.0:
        raise StructureError(f"plddt[{index}] must be finite and in 0..100")
    return number


def _pae(rows: list[Any], length: int) -> tuple[tuple[float, ...], ...]:
    if len(rows) != length:
        raise StructureError("pae must be square and match the sequence")
    parsed: list[tuple[float, ...]] = []
    for row in rows:
        if not isinstance(row, list) or len(row) != length:
            raise StructureError("pae must be square and match the sequence")
        parsed.append(tuple(_angstrom(value) for value in row))
    return tuple(parsed)


def _angstrom(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise StructureError("pae must contain only finite numbers")
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise StructureError("pae must contain only finite, non-negative numbers")
    return number


def _mean_pae(matrix: tuple[tuple[float, ...], ...]) -> float:
    total = sum(sum(row) for row in matrix)
    count = sum(len(row) for row in matrix)
    return total / count
