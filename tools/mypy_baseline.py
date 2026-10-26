#!/usr/bin/env python3
"""Fail when a module has more mypy errors than the baseline allows.

``tools/mypy-baseline.txt`` lists modules that already have errors, one per
line: a path relative to the repository root, then the allowed count.
A module that is not listed must be clean. A listed module may not grow
errors. Counts that have fallen are reported and fail too, so the baseline
stays the current number rather than a ceiling nobody lowers.

Usage, from the repository root::

    python tools/mypy_baseline.py

The package must be installed (``pip install -e ".[typecheck]"``). Only the
standard library and mypy are used.
"""

from __future__ import annotations

import contextlib
import subprocess
import sys
from collections import Counter
from pathlib import Path

BASELINE = Path("tools/mypy-baseline.txt")


def parse_baseline(text: str) -> dict[str, int]:
    """Path to allowed error count. Blank lines and ``#`` comments are ignored."""
    allowed: dict[str, int] = {}
    for line_number, raw in enumerate(text.splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        path, _, count = line.partition(" ")
        if not count or not count.isdigit() or " " in count:
            raise ValueError(f"{BASELINE}:{line_number}: expected 'path count', got {raw!r}")
        allowed[path] = int(count)
    return allowed


def errors_by_module(report: str, root: Path) -> Counter[str]:
    """Count mypy ``error:`` lines per file, as paths relative to ``root``."""
    counts: Counter[str] = Counter()
    for line in report.splitlines():
        if ": error:" not in line:
            continue
        path = line.split(":", 1)[0].strip()
        if not path:
            continue
        as_path = Path(path)
        if as_path.is_absolute():
            with contextlib.suppress(ValueError):
                as_path = as_path.relative_to(root)
        counts[as_path.as_posix()] += 1
    return counts


def compare(actual: Counter[str], allowed: dict[str, int]) -> list[str]:
    """Human-readable failures. Empty when the tree matches the baseline."""
    problems: list[str] = []
    for path in sorted(set(actual) | set(allowed)):
        have = actual.get(path, 0)
        if path not in allowed:
            if have:
                problems.append(f"{path}: {have} error(s), and it is not in the baseline")
            continue
        permitted = allowed[path]
        if have > permitted:
            problems.append(f"{path}: {have} error(s), baseline allows {permitted}")
        elif have < permitted:
            problems.append(f"{path}: {have} error(s), baseline still allows {permitted}; lower it")
    return problems


def run_mypy(root: Path) -> str:
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "mypy",
            "--package",
            "marse",
            "--no-error-summary",
            "--show-error-codes",
        ],
        cwd=root,
        capture_output=True,
        text=True,
    )
    # mypy exits 1 when it found errors. That is the input to the baseline, not
    # a failure of the tool. Other exits mean mypy itself did not run.
    if completed.returncode not in (0, 1):
        sys.stderr.write(completed.stderr)
        raise SystemExit(completed.returncode or 1)
    return completed.stdout


def main(argv: list[str] | None = None) -> int:
    root = Path(argv[0]).resolve() if argv else Path.cwd()
    baseline_path = root / "tools" / "mypy-baseline.txt"
    allowed = parse_baseline(baseline_path.read_text(encoding="utf-8"))
    actual = errors_by_module(run_mypy(root), root)
    problems = compare(actual, allowed)
    if problems:
        print("\n".join(problems))
        print(f"\nmypy baseline: {len(problems)} module(s) differ from {baseline_path}")
        return 1
    print(f"mypy baseline: {sum(actual.values())} known error(s) in {len(allowed)} module(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
