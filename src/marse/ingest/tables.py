"""Read a CSV or TSV table. Values stay text until an importer checks them."""

from __future__ import annotations

import csv
import io
from pathlib import Path

from marse.ingest.errors import IngestError


def read_table(path: Path) -> tuple[tuple[str, ...], tuple[dict[str, str], ...]]:
    """Header and rows. A tab-separated file is read as TSV; anything else as CSV."""
    text = path.read_text(encoding="utf-8")
    if not text.strip():
        raise IngestError(f"{path.name}: empty table")
    first = text.splitlines()[0]
    delimiter = "\t" if "\t" in first else ","
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter, restkey="_extra")
    if reader.fieldnames is None:
        raise IngestError(f"{path.name}: a table needs a header")
    columns = tuple(name.strip() for name in reader.fieldnames if name is not None)
    if any(not name for name in columns) or len(columns) != len(set(columns)):
        raise IngestError(f"{path.name}: column names must be present and unique")
    rows: list[dict[str, str]] = []
    for index, raw in enumerate(reader, start=2):
        if raw.get("_extra"):
            raise IngestError(f"{path.name}:{index}: more fields than the header")
        row = {key.strip(): (value or "").strip() for key, value in raw.items() if key != "_extra"}
        if not any(row.values()):
            continue
        rows.append(row)
    if not rows:
        raise IngestError(f"{path.name}: no data rows")
    return columns, tuple(rows)


def require_columns(path_name: str, columns: tuple[str, ...], required: tuple[str, ...]) -> None:
    missing = [name for name in required if name not in columns]
    if missing:
        raise IngestError(f"{path_name}: missing column {missing[0]}")
