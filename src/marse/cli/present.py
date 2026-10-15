"""Print the lines ``marse run`` and ``marse ecosystem`` have always printed."""

from __future__ import annotations

from marse.cli.summarise import _summarise, _summarise_in_space, _summarise_well_mixed
from marse.core.reactive_transport import ReactiveTransportResult
from marse.core.simulation import BiofilmProfileResult, SimulationResult
from marse.core.well_mixed import WellMixedResult
from marse.ecosystem import EcosystemResult
from marse.engines.base import Completed
from marse.engines.legacy_ecosystem import UNVERIFIED_WARNING

__all__ = [
    "UNVERIFIED_WARNING",
    "present_batch",
    "present_ecosystem",
    "present_reactive",
    "present_well_mixed",
    "warn_unverified",
]


def warn_unverified() -> None:
    """Tell the operator that ecosystem v1 results are not a validated model."""
    import sys

    print(UNVERIFIED_WARNING, file=sys.stderr)


def present_batch(completed: Completed) -> None:
    result = completed.result
    if not isinstance(result, SimulationResult | BiofilmProfileResult):
        raise TypeError(f"batch summary expected a batch result, got {type(result).__name__}")
    _summarise(result)
    _print_pair(completed.written)


def present_well_mixed(completed: Completed) -> None:
    result = completed.result
    if not isinstance(result, WellMixedResult):
        raise TypeError(
            f"well-mixed summary expected a well-mixed result, got {type(result).__name__}"
        )
    _summarise_well_mixed(result)
    _print_pair(completed.written)


def present_reactive(completed: Completed) -> None:
    result = completed.result
    if not isinstance(result, ReactiveTransportResult):
        raise TypeError(
            f"spatial summary expected a reactive-transport result, got {type(result).__name__}"
        )
    _summarise_in_space(result)
    print()
    for path in completed.written:
        print(f"wrote {path}")
    pvd = next(path for path in completed.written if path.suffix == ".pvd")
    print(f"\nopen {pvd} in ParaView to explore the run in 3-D")
    print(f"replay it with:  marse replay {completed.written[-1]}")


def present_ecosystem(completed: Completed) -> None:
    result = completed.result
    if not isinstance(result, EcosystemResult):
        raise TypeError(
            f"ecosystem summary expected an ecosystem result, got {type(result).__name__}"
        )
    manifest = result.manifest
    if manifest is None:
        raise RuntimeError("the ecosystem engine returned no manifest")
    frames, viewer, manifest_path = completed.written
    print(f"run_id      {manifest.run_id}")
    print(f"experiment  {result.config.experiment_id}")
    print(f"steps       {result.final_state.step} over {result.final_state.time_h:g} h")
    print(
        f"frames      {completed.stored_frames} stored, one every {completed.frame_every} step(s)"
    )
    print(f"wrote       {frames}")
    print(f"wrote       {viewer}")
    print(f"wrote       {manifest_path}")
    print(f"\nreplay it with:  marse replay {manifest_path}")


def _print_pair(written: tuple) -> None:
    trajectory, manifest = written
    print(f"\nwrote {trajectory}")
    print(f"wrote {manifest}")
    print(f"\nreplay it with:  marse replay {manifest}")
