"""``marse bench``: list cases, run them, and print a scored report.

``list`` does not run. ``run`` writes ``bench.json`` and ``bench.html`` beside
the run directory. ``report`` reads those files and does not re-run.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from marse.bench.adapters import adapter_for
from marse.bench.cases import Case, CaseError, load_cases
from marse.bench.metrics import evaluate
from marse.bench.report import document, read_reports, write_report
from marse.cli.present import warn_unverified


def command_bench(args: argparse.Namespace) -> int:
    """``marse bench`` with no subcommand."""
    print("marse: bench needs list, run, or report")
    return 2


def command_bench_list(args: argparse.Namespace) -> int:
    cases = load_cases(Path(args.cases), Path(args.root).resolve())
    if not cases:
        print(f"no cases in {args.cases}")
        return 0
    for case in cases:
        print(f"{case.id}  {case.adapter}  {case.experiment}")
    return 0


def command_bench_run(args: argparse.Namespace) -> int:
    root = Path(args.root).resolve()
    cases = load_cases(Path(args.cases), root)
    selected = _select(cases, args.case)
    if not selected:
        raise CaseError(f"no cases in {args.cases}")
    output = Path(args.output)
    failed = False
    for case in selected:
        adapter = adapter_for(case.adapter)
        case_dir = output / case.id
        record = adapter.run(case, case_dir, root)
        if record.unverified:
            warn_unverified()
        evaluation = evaluate(case, case_dir, unverified=record.unverified)
        payload = document(
            case,
            evaluation,
            engine=record.engine,
            unverified=record.unverified,
            container=record.container,
        )
        write_report(case_dir, payload)
        status = "pass" if payload["passed"] else "fail"
        print(f"{case.id}  {status}")
        failed = failed or not payload["passed"]
    return 1 if failed else 0


def command_bench_report(args: argparse.Namespace) -> int:
    rows = read_reports(Path(args.directory))
    if not rows:
        raise CaseError(f"no bench.json under {args.directory}")
    failed = False
    for row in rows:
        status = "pass" if row["passed"] else "fail"
        print(f"{row['case_id']}  {status}  {row['engine']}")
        if row.get("reason"):
            print(f"  {row['reason']}")
        failed = failed or not row["passed"]
    return 1 if failed else 0


def _select(cases: tuple[Case, ...], wanted: list[str] | None) -> tuple[Case, ...]:
    if not wanted:
        return cases
    by_id = {case.id: case for case in cases}
    missing = [name for name in wanted if name not in by_id]
    if missing:
        raise CaseError(f"unknown case: {missing[0]}")
    return tuple(by_id[name] for name in wanted)
