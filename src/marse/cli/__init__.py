"""Command-line entry point: ``marse``.

The commands are:

    marse run experiment.json        run a batch, biofilm, or version-2 network
    marse ecosystem experiment.json  run a 2-D multispecies ecosystem (unverified)
    marse replay manifest.json       re-run any of them from its manifest and compare
    marse ls runs/                   list run directories
    marse compare runs/a runs/b      compare two manifests
    marse check network.json         check a reaction network (schema version 2)
    marse bench list cases/          list benchmark cases without running them
    marse bench run cases/ -o bench/ run the cases and write a scored report
    marse bench report bench/        print a scored report without re-running

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
from marse.bench.cases import CaseError
from marse.cli.bench import (
    command_bench,
    command_bench_list,
    command_bench_report,
    command_bench_run,
)
from marse.cli.check import _command_check
from marse.cli.ecosystem import command_ecosystem, command_ecosystem_batch
from marse.cli.niche import command_niche_scan
from marse.cli.replay import command_replay
from marse.cli.run import command_run
from marse.cli.runs import command_compare, command_ls
from marse.core.config import ConfigError
from marse.engines.legacy_ecosystem import MAX_STORED_FRAMES
from marse.processes.niche import NicheError

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

    ls_command = commands.add_parser("ls", help="list runs under a directory")
    ls_command.add_argument("directory", help="a run directory, or a directory that contains runs")
    ls_command.add_argument("--json", action="store_true", help="print the index as JSON")
    ls_command.set_defaults(handler=command_ls)

    compare_command = commands.add_parser(
        "compare", help="compare two manifests without re-running them"
    )
    compare_command.add_argument("left", help="a manifest.json, or a run directory")
    compare_command.add_argument("right", help="a manifest.json, or a run directory")
    compare_command.set_defaults(handler=command_compare)

    bench_command = commands.add_parser("bench", help="list, run, or report benchmark cases")
    bench_command.set_defaults(handler=command_bench)
    bench_sub = bench_command.add_subparsers(dest="bench_command")

    bench_list = bench_sub.add_parser("list", help="list cases without running them")
    bench_list.add_argument("cases", help="directory of YAML case files")
    bench_list.add_argument("--root", default=".", help="repository root for experiment paths")
    bench_list.set_defaults(handler=command_bench_list)

    bench_run = bench_sub.add_parser("run", help="run cases and write a scored report")
    bench_run.add_argument("cases", help="directory of YAML case files")
    bench_run.add_argument(
        "-o", "--output", required=True, help="directory for run and report files"
    )
    bench_run.add_argument("--root", default=".", help="repository root for experiment paths")
    bench_run.add_argument(
        "--case", action="append", help="run only this case id; repeat to select more"
    )
    bench_run.set_defaults(handler=command_bench_run)

    bench_report = bench_sub.add_parser("report", help="print a scored report without re-running")
    bench_report.add_argument("directory", help="directory written by 'marse bench run'")
    bench_report.set_defaults(handler=command_bench_report)
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
    except CaseError as error:
        print(f"marse: invalid benchmark case: {error}")
        return 2
    except (OSError, ValueError, KeyError) as error:
        print(f"marse: {type(error).__name__}: {error}")
        return 2
    finally:
        np.seterr(**errors)  # NumPy's error handling as the caller had it
