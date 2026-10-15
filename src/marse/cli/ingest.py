"""``marse ingest``: check a table and write a feed. It does not run a model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from marse.ingest import importer_for, registered_kinds


def command_ingest(args: argparse.Namespace) -> int:
    feed = importer_for(args.kind).load(Path(args.table))
    text = json.dumps(feed.to_dict(), indent=2, sort_keys=True) + "\n"
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")
    print(f"{feed.kind}  {len(feed.rows)} rows  qc pass")
    return 0


def add_ingest_parser(commands: Any) -> None:
    """Register ``marse ingest`` on the root parser."""
    ingest = commands.add_parser("ingest", help="check a table and write a feed")
    ingest.add_argument(
        "kind",
        choices=registered_kinds(),
        help="abundance, transcriptomics, metabolomics, or environment",
    )
    ingest.add_argument("table", help="CSV or TSV table")
    ingest.add_argument("-o", "--output", required=True, help="path of the feed JSON")
    ingest.set_defaults(handler=command_ingest)
