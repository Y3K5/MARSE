"""``marse ecosystem`` and ``marse ecosystem-batch``.

The two-dimensional engine is unverified. Every invocation warns on stderr
before it runs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from marse.analysis.ensemble import ScenarioBatch, run_batch
from marse.cli.present import present_ecosystem, warn_unverified
from marse.engines import registry


def command_ecosystem(args: argparse.Namespace) -> int:
    warn_unverified()
    engine = registry.get("legacy_ecosystem_v1")
    source = Path(args.experiment)
    plan = engine.load(source)
    output_dir = Path(args.output) if args.output else engine.default_output_dir(plan, source)
    completed = engine.run(plan, output_dir, frame_interval_h=args.frame_interval_h)
    present_ecosystem(completed)
    return 0


def command_ecosystem_batch(args: argparse.Namespace) -> int:
    warn_unverified()
    raw = json.loads(Path(args.experiment).read_text(encoding="utf-8"))
    batch = ScenarioBatch.grid(raw["base"], raw.get("parameters", {}))
    catalog = run_batch(batch, workers=args.workers)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    catalog_path = catalog.write_json(output / "catalog.json")
    print(f"scenarios   {len(batch.scenarios)}")
    print(f"completed   {len(catalog.query(status='completed'))}")
    print(f"failed      {len(catalog.query(status='failed'))}")
    print(f"wrote       {catalog_path}")
    return 0
