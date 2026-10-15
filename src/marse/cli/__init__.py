"""Command-line entry point: ``marse``.

The commands are unchanged:

    marse run experiment.json        run a batch, biofilm, or version-2 network
    marse ecosystem experiment.json  run a 2-D multispecies ecosystem (unverified)
    marse replay manifest.json       re-run any of them from its manifest and compare
    marse check network.json         check a reaction network (schema version 2)

``run`` selects the engine from the file. A top-level ``engine`` field names
one explicitly (``legacy_batch``, ``steady_biofilm``, ``well_mixed``,
``reactive_transport``, ``legacy_ecosystem_v1``). Without it, the choice is
the one these commands have always made. ``replay`` is the reproducibility
claim made executable: it rebuilds the configuration from the manifest,
verifies the checksum, runs again, and reports whether the results match.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence

import numpy as np

from marse import __version__
from marse.cli.check import _command_check
from marse.cli.ecosystem import command_ecosystem, command_ecosystem_batch
from marse.cli.niche import command_niche_scan
from marse.cli.replay import command_replay
from marse.cli.run import command_run
from marse.core.config import ConfigError
from marse.engines.legacy_ecosystem import MAX_STORED_FRAMES
from marse.microbes.niche import NicheError

__all__ = ["build_parser", "main"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="marse",
        description="MARSE: Microbial Adaptability Resource Simulation Engine.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command")

    run_command = commands.add_parser("run", help="run a simulation from an experiment file")
    run_command.add_argument("experiment", help="path to a JSON experiment configuration")
    run_command.add_argument(
        "-o", "--output", help="directory for outputs (default: runs/<run_id>)"
    )
    run_command.set_defaults(handler=command_run)

    replay_command = commands.add_parser("replay", help="re-run from a manifest and compare")
    replay_command.add_argument("manifest", help="path to a manifest.json written by 'marse run'")
    replay_command.add_argument("-o", "--output", help="directory for the replayed outputs")
    replay_command.set_defaults(handler=command_replay)

    ecosystem_command = commands.add_parser(
        "ecosystem",
        help="run a 2D multi-species ecosystem and export a browser viewer (unverified)",
    )
    ecosystem_command.add_argument("experiment", help="path to an ecosystem JSON configuration")
    ecosystem_command.add_argument("-o", "--output", help="directory for frames and viewer")
    ecosystem_command.add_argument(
        "--frame-interval-h",
        type=float,
        help=f"hours between stored frames (default: about {MAX_STORED_FRAMES} frames per run)",
    )
    ecosystem_command.set_defaults(handler=command_ecosystem)

    niche_command = commands.add_parser(
        "niche-scan", help="scan environmental conditions against species capabilities"
    )
    niche_command.add_argument("experiment", help="path to a niche scan JSON configuration")
    niche_command.add_argument("-o", "--output", help="directory for scan outputs")
    niche_command.set_defaults(handler=command_niche_scan)

    batch_command = commands.add_parser(
        "ecosystem-batch", help="run a reproducible ecosystem parameter ensemble"
    )
    batch_command.add_argument("experiment", help="JSON file with base and parameters")
    batch_command.add_argument("-o", "--output", required=True, help="directory for catalog output")
    batch_command.add_argument("--workers", type=int, default=1, help="parallel worker processes")
    batch_command.set_defaults(handler=command_ecosystem_batch)

    check_command = commands.add_parser(
        "check", help="check a reaction network and print its balanced processes"
    )
    check_command.add_argument("network", help="path to a network JSON file (schema version 2)")
    check_command.set_defaults(handler=_command_check)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return 0
    errors = np.geterr()
    try:
        return int(args.handler(args))
    except ConfigError as error:
        print(f"marse: invalid experiment: {error}")
        return 2
    except NicheError as error:
        print(f"marse: invalid niche scan: {error}")
        return 2
    except (OSError, ValueError, KeyError) as error:
        print(f"marse: {type(error).__name__}: {error}")
        return 2
    finally:
        np.seterr(**errors)  # NumPy's error handling as the caller had it
