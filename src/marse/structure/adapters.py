"""ESMFold and AlphaFold adapters. Each calls a function it was given.

Neither adapter imports a folding library. A deployment that has one passes
the function that runs it. A missing function is refused before any score is
read.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any, Protocol

from marse.structure.record import StructureError, StructureRecord, score_prediction

Payload = Mapping[str, Any]
Caller = Callable[[str], Payload]


class FoldAdapter(Protocol):
    """Turn one sequence into one checked structure record."""

    name: str

    def predict(self, sequence: str) -> StructureRecord:
        """Return a record, or raise :class:`StructureError`."""


class ServiceAdapter:
    """An adapter whose network or local model call is injected."""

    def __init__(self, name: str, caller: Caller | None) -> None:
        self.name = name
        self._caller = caller

    def predict(self, sequence: str) -> StructureRecord:
        if self._caller is None:
            raise StructureError(f"{self.name} caller is not configured")
        payload = self._caller(sequence)
        if not isinstance(payload, Mapping):
            raise StructureError(f"{self.name} returned no score mapping")
        plddt = payload.get("plddt")
        if not isinstance(plddt, list):
            raise StructureError(f"{self.name} payload needs a plddt list")
        pae = payload.get("pae")
        if pae is not None and not isinstance(pae, list):
            raise StructureError(f"{self.name} pae must be a matrix when it is present")
        return score_prediction(self.name, sequence, plddt, pae)


def esmfold_adapter(caller: Caller | None = None) -> ServiceAdapter:
    """ESMFold scores. pLDDT is required. PAE may be absent."""
    return ServiceAdapter("esmfold", caller)


def alphafold_adapter(caller: Caller | None = None) -> ServiceAdapter:
    """AlphaFold scores. pLDDT is required. PAE is recorded when the caller sends it."""
    return ServiceAdapter("alphafold", caller)
