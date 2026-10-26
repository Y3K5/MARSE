"""Write the data file and manifest a run leaves behind.

The data file depends on the result: a depth profile, spatial totals, or a
well-mixed trajectory. The manifest is always ``manifest.json``.
"""

from __future__ import annotations

from pathlib import Path

from marse.core.reactive_transport import ReactiveTransportResult
from marse.core.simulation import BiofilmProfileResult, SimulationResult
from marse.core.well_mixed import WellMixedResult
from marse.io.layout import publish_run

Result = SimulationResult | BiofilmProfileResult | WellMixedResult | ReactiveTransportResult


def write_run_files(result: Result, output_dir: Path) -> tuple[Path, Path]:
    """Write the run's data file and its manifest. The data file differs by kind."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if isinstance(result, BiofilmProfileResult):
        data = result.write_profile(output_dir / "profile.csv")
    elif isinstance(result, ReactiveTransportResult):
        data = result.write_totals(output_dir / "totals.csv")
    else:
        data = result.write_trajectory(output_dir / "trajectory.csv")
    manifest = result.manifest.write(output_dir / "manifest.json")
    publish_run(output_dir, result.manifest)
    return data, manifest
