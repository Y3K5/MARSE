"""Adapters that produce one run directory.

The MARSE adapter calls the engine registry in-process. Its container field is
empty: it does not start a shell and it does not start a container. A later
adapter may set that field to an image name without changing the scorer.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from marse.bench.cases import Case, CaseError
from marse.engines import registry
from marse.engines.legacy_ecosystem import legacy_ecosystem_v1


@dataclass(frozen=True, slots=True)
class RunRecord:
    """The three files every engine already writes, plus who wrote them."""

    manifest_path: Path
    summary_path: Path
    report_path: Path
    engine: str
    unverified: bool
    container: str


class MarseAdapter:
    """Run a case through the engine registry, in this process."""

    name = "marse"
    container = ""

    def run(self, case: Case, output_dir: Path, root: Path) -> RunRecord:
        experiment = root / case.experiment
        output_dir.mkdir(parents=True, exist_ok=True)
        if case.command == "none":
            return RunRecord(
                output_dir / "manifest.json",
                output_dir / "summary.json",
                output_dir / "report.html",
                "not_run",
                False,
                self.container,
            )
        if case.command == "ecosystem":
            plan = legacy_ecosystem_v1.load(experiment)
            legacy_ecosystem_v1.run(plan, output_dir)
            name = legacy_ecosystem_v1.name
            unverified = legacy_ecosystem_v1.unverified
        else:
            engine, plan = registry.load(experiment)
            engine.run(plan, output_dir)
            name = engine.name
            unverified = engine.unverified
        return RunRecord(
            output_dir / "manifest.json",
            output_dir / "summary.json",
            output_dir / "report.html",
            name,
            bool(unverified),
            self.container,
        )


def adapter_for(name: str) -> MarseAdapter:
    """The adapter a case names. Only the in-process MARSE adapter is registered."""
    if name != MarseAdapter.name:
        raise CaseError(f"unknown adapter: {name}")
    return MarseAdapter()
