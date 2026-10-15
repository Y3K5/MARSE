"""One importer per table kind. Registration is the extension point."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from marse.ingest.errors import IngestError
from marse.ingest.feed import Feed


class Importer(Protocol):
    """Read one kind of table into a :class:`~marse.ingest.feed.Feed`."""

    @property
    def kind(self) -> str:
        """The table kind this importer accepts."""

    def load(self, path: Path) -> Feed:
        """Check ``path`` and return the feed. Raise :class:`IngestError` when it fails."""


_REGISTRY: dict[str, Importer] = {}


def register(importer: Importer) -> None:
    """Add an importer. A second registration of the same kind is refused."""
    if importer.kind in _REGISTRY:
        raise IngestError(f"importer already registered: {importer.kind}")
    _REGISTRY[importer.kind] = importer


def unregister(kind: str) -> None:
    """Drop one importer. Tests use this so a stand-in does not stay registered."""
    _REGISTRY.pop(kind, None)


def importer_for(kind: str) -> Importer:
    """The importer registered for ``kind``."""
    try:
        return _REGISTRY[kind]
    except KeyError:
        raise IngestError(f"unknown importer: {kind}") from None


def registered_kinds() -> tuple[str, ...]:
    """Kinds that can be loaded, in alphabetical order."""
    return tuple(sorted(_REGISTRY))


def _number(text: str, where: str) -> float:
    try:
        value = float(text)
    except ValueError:
        raise IngestError(f"{where}: expected a number") from None
    if not math.isfinite(value):
        raise IngestError(f"{where}: must be finite")
    return value


_RELATIVE_SUM = 1e-6


@dataclass(frozen=True, slots=True)
class EnvironmentImporter:
    """``name,value,unit`` rows that match the batch environment fields."""

    kind: str = "environment"

    def load(self, path: Path) -> Feed:
        from marse.ingest.tables import read_table, require_columns

        columns, rows = read_table(path)
        require_columns(path.name, columns, ("name", "value", "unit"))
        environment: dict[str, float] = {}
        allowed = {"temperature_c": "C", "ph": "pH"}
        for index, row in enumerate(rows, start=2):
            name = row["name"]
            if name not in allowed:
                raise IngestError(f"{path.name}:{index}: {name} is not an environment field")
            if row["unit"] != allowed[name]:
                raise IngestError(f"{path.name}:{index}: {name} unit must be {allowed[name]}")
            if name in environment:
                raise IngestError(f"{path.name}:{index}: {name} is repeated")
            value = _number(row["value"], f"{path.name}:{index}")
            if name == "temperature_c" and not -20.0 <= value <= 150.0:
                raise IngestError(f"{path.name}:{index}: temperature_c outside -20..150 C")
            if name == "ph" and not 0.0 <= value <= 14.0:
                raise IngestError(f"{path.name}:{index}: ph outside 0..14")
            environment[name] = value
        if set(environment) != set(allowed):
            raise IngestError(f"{path.name}: temperature_c and ph are both required")
        return Feed(
            self.kind,
            path.name,
            columns,
            rows,
            {"environment": environment},
            "temperature_c and ph are copied. No other environment field is.",
        )


@dataclass(frozen=True, slots=True)
class AbundanceImporter:
    """``sample,taxon,abundance,unit`` with relative fractions or raw counts."""

    kind: str = "abundance"

    def load(self, path: Path) -> Feed:
        from marse.ingest.tables import read_table, require_columns

        columns, rows = read_table(path)
        require_columns(path.name, columns, ("sample", "taxon", "abundance", "unit"))
        seen: set[tuple[str, str]] = set()
        totals: dict[str, float] = {}
        for index, row in enumerate(rows, start=2):
            if not row["sample"] or not row["taxon"]:
                raise IngestError(f"{path.name}:{index}: sample and taxon are required")
            key = (row["sample"], row["taxon"])
            if key in seen:
                raise IngestError(
                    f"{path.name}:{index}: {row['taxon']} is repeated in {row['sample']}"
                )
            seen.add(key)
            if row["unit"] not in {"relative", "count"}:
                raise IngestError(f"{path.name}:{index}: unit must be relative or count")
            value = _number(row["abundance"], f"{path.name}:{index}")
            if value < 0:
                raise IngestError(f"{path.name}:{index}: abundance must be non-negative")
            if row["unit"] == "relative":
                totals[row["sample"]] = totals.get(row["sample"], 0.0) + value
        for sample, total in totals.items():
            if abs(total - 1.0) > _RELATIVE_SUM:
                raise IngestError(f"{path.name}: relative abundances in {sample} sum to {total}")
        return Feed(
            self.kind,
            path.name,
            columns,
            rows,
            None,
            "Abundances are not biomass and are not written into an experiment.",
        )


@dataclass(frozen=True, slots=True)
class _FeatureImporter:
    """A sample/feature/value/unit table. Transcripts never become a rate."""

    kind: str
    units: frozenset[str]
    limits: str
    fragment: bool

    def load(self, path: Path) -> Feed:
        from marse.ingest.tables import read_table, require_columns

        columns, rows = read_table(path)
        require_columns(path.name, columns, ("sample", "feature", "value", "unit"))
        seen: set[tuple[str, str]] = set()
        parsed: list[tuple[str, str, float, str]] = []
        for index, row in enumerate(rows, start=2):
            if not row["sample"] or not row["feature"]:
                raise IngestError(f"{path.name}:{index}: sample and feature are required")
            key = (row["sample"], row["feature"])
            if key in seen:
                raise IngestError(
                    f"{path.name}:{index}: {row['feature']} is repeated in {row['sample']}"
                )
            seen.add(key)
            if row["unit"] not in self.units:
                raise IngestError(f"{path.name}:{index}: unit {row['unit']!r} is not accepted")
            value = _number(row["value"], f"{path.name}:{index}")
            if value < 0:
                raise IngestError(f"{path.name}:{index}: value must be non-negative")
            parsed.append((row["sample"], row["feature"], value, row["unit"]))
        return Feed(
            self.kind,
            path.name,
            columns,
            rows,
            _substrate_fragment(parsed) if self.fragment else None,
            self.limits,
        )


def _substrate_fragment(
    rows: list[tuple[str, str, float, str]],
) -> dict[str, dict[str, float | str]] | None:
    """One metabolite, one sample, in mol per m3, becomes one substrate initial amount.

    1 mmol/L is 1 mol/m3, which is the batch field ``initial_mm``. Two features,
    two samples, or another unit are stored and not copied.
    """
    samples = {row[0] for row in rows}
    features = {row[1] for row in rows}
    units = {row[3] for row in rows}
    if len(samples) != 1 or len(features) != 1 or len(rows) != 1:
        return None
    if units != {"mmol_per_l"} and units != {"mol_per_m3"}:
        return None
    _sample, feature, value, _unit = rows[0]
    return {"substrate": {"name": feature, "initial_mm": value}}


def _register_builtins() -> None:
    register(EnvironmentImporter())
    register(AbundanceImporter())
    register(
        _FeatureImporter(
            "transcriptomics",
            frozenset({"count", "tpm"}),
            "Counts and TPM are not rates and are not written into an experiment.",
            False,
        )
    )
    register(
        _FeatureImporter(
            "metabolomics",
            frozenset({"mmol_per_l", "mol_per_m3", "arbitrary"}),
            "One metabolite in mol per m3 can fill substrate.initial_mm. "
            "That is the table's number, not a yield or a calibrated uptake.",
            True,
        )
    )


_register_builtins()
