"""``marse ls`` and ``marse compare``: the run index, read back.

``ls`` lists manifests under a directory. ``compare`` reads two manifests and
reports whether the configuration checksum and the recorded outputs match.
Neither command re-runs a simulation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from marse.io.layout import list_runs, read_manifest


def command_ls(args: argparse.Namespace) -> int:
    rows = list_runs(Path(args.directory))
    if not rows:
        print(f"no runs under {args.directory}")
        return 0
    if args.json:
        print(json.dumps(rows, indent=2, sort_keys=True))
        return 0
    print(f"{'directory':<28} {'kind':<20} {'experiment':<28} checksum")
    for row in rows:
        print(
            f"{row['directory']:<28} {row['kind']:<20} {row['experiment_id']:<28} "
            f"{row['config_sha256'][:12]}"
        )
    return 0


def command_compare(args: argparse.Namespace) -> int:
    left = read_manifest(Path(args.left))
    right = read_manifest(Path(args.right))
    same_config = left.config_sha256 == right.config_sha256
    same_outputs = left.outputs == right.outputs
    print(f"left        {left.experiment_id}  {left.kind}  {left.config_sha256[:16]}")
    print(f"right       {right.experiment_id}  {right.kind}  {right.config_sha256[:16]}")
    print(f"config      {'same checksum' if same_config else 'DIFFERENT checksum'}")
    outputs = "same recorded outputs" if same_outputs else "DIFFERENT recorded outputs"
    print(f"outputs     {outputs}")
    if left.kind != right.kind:
        print(f"kind        {left.kind} vs {right.kind}")
    return 0 if same_config and same_outputs else 1
