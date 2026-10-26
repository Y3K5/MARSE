"""``marse niche-scan``: environmental conditions against species capabilities."""

from __future__ import annotations

import argparse
from pathlib import Path

from marse.processes.niche import load_niche_scan, run_niche_scan


def command_niche_scan(args: argparse.Namespace) -> int:
    scan = load_niche_scan(args.experiment)
    result = run_niche_scan(scan)
    output_dir = (
        Path(args.output)
        if args.output
        else Path(args.experiment).parent / "runs" / scan.experiment_id
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = result.write_json(output_dir / "niche-scan.json")
    csv_path = result.write_csv(output_dir / "niche-scan.csv")
    print(f"experiment  {scan.experiment_id}")
    print(f"rows        {len(result.rows)}")
    print(f"wrote       {json_path}")
    print(f"wrote       {csv_path}")
    return 0
