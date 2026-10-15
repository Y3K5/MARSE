"""``marse run``: choose an engine from the experiment file and run it."""

from __future__ import annotations

import argparse
from pathlib import Path

from marse.cli.present import (
    present_batch,
    present_ecosystem,
    present_reactive,
    present_well_mixed,
    warn_unverified,
)
from marse.engines import registry

_PRESENT = {
    "legacy_batch": present_batch,
    "steady_biofilm": present_batch,
    "well_mixed": present_well_mixed,
    "reactive_transport": present_reactive,
    "legacy_ecosystem_v1": present_ecosystem,
}


def command_run(args: argparse.Namespace) -> int:
    source = Path(args.experiment)
    engine, plan = registry.load(source)
    if engine.unverified:
        warn_unverified()
    output_dir = Path(args.output) if args.output else engine.default_output_dir(plan, source)
    completed = engine.run(plan, output_dir)
    _PRESENT[engine.name](completed)
    return 0
