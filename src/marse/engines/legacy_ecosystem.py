"""The version-1 two-dimensional ecosystem engine.

Unverified. ``docs/validation.md`` records the known defects, and the command
line warns before every run and replay. The version-1 periodontal file still
uses this engine. A separate version-2 column lives in
``examples/networks/periodontal_pathogen_v2.json`` and does not run here.
The implementation remains :mod:`marse.ecosystem`.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from marse.core.provenance import Manifest
from marse.ecosystem import EcosystemConfig, EcosystemResult, write_viewer
from marse.ecosystem import load_experiment as load_ecosystem_experiment
from marse.ecosystem import run as run_ecosystem
from marse.ecosystem.framestore import EcosystemFrameSink, FrameStore, StoredFrames
from marse.ecosystem.model import ecosystem_from_dict, recorded_steps
from marse.engines.base import Completed
from marse.io.layout import publish_run

MAX_STORED_FRAMES = 200
"""About how many frames ``marse ecosystem`` stores by default, whatever the run length."""

UNVERIFIED_WARNING = (
    "marse: warning: the 2-D ecosystem engine (ecosystem v1) is unverified. "
    "docs/validation.md records known defects (mass balance, oxygen transport, "
    "and others). Results, including the periodontal study, support no "
    "quantitative conclusion."
)


class LegacyEcosystem:
    """Cell-based two-dimensional ecosystem. Not a validated model."""

    name = "legacy_ecosystem_v1"
    kinds: tuple[str, ...] = ("ecosystem",)
    unverified = True
    writes_during_replay = True

    def load(self, path: Path) -> EcosystemConfig:
        return load_ecosystem_experiment(path)

    def validate(self, document: dict[str, Any]) -> EcosystemConfig:
        return ecosystem_from_dict(document)

    def default_output_dir(self, plan: Any, source: Path) -> Path:
        if not isinstance(plan, EcosystemConfig):
            raise TypeError("legacy_ecosystem_v1 runs an EcosystemConfig")
        return source.parent / "runs" / plan.experiment_id

    def run(
        self,
        plan: Any,
        output_dir: Path,
        *,
        frame_interval_h: float | None = None,
    ) -> Completed:
        if not isinstance(plan, EcosystemConfig):
            raise TypeError("legacy_ecosystem_v1 runs an EcosystemConfig")
        frame_every = _frame_every(plan, frame_interval_h)
        result, stored, written = _run_into(plan, output_dir, frame_every)
        return Completed(result, written, stored_frames=stored, frame_every=frame_every)

    def replay(self, manifest: Manifest, output_dir: Path | None) -> Completed:
        config = manifest.experiment()
        if not isinstance(config, EcosystemConfig):
            raise ValueError(
                f"manifest declares {manifest.kind!r} but its configuration is not an ecosystem"
            )
        if output_dir is None:
            return Completed(run_ecosystem(config, frame_every=None))
        frame_every = _frame_every(config, None)
        result, stored, written = _run_into(config, output_dir, frame_every)
        return Completed(result, written, stored_frames=stored, frame_every=frame_every)

    def write_outputs(self, result: Any, output_dir: Path) -> tuple[Path, ...]:
        raise RuntimeError("an ecosystem replay writes its files while it runs")


def _frame_every(config: EcosystemConfig, interval_h: float | None) -> int:
    """Steps between stored frames: from ``--frame-interval-h``, else a cap on the count."""
    if interval_h is not None:
        if not interval_h > 0:
            raise ValueError("--frame-interval-h must be positive")
        return max(1, round(interval_h / config.timestep_h))
    return max(1, math.ceil(config.steps / MAX_STORED_FRAMES))


def _run_into(
    config: EcosystemConfig, output_dir: Path, frame_every: int
) -> tuple[EcosystemResult, int, tuple[Path, ...]]:
    """Run, streaming frames into ``output_dir/frames``; then write the viewer and manifest.

    Memory stays flat however long the run: frames go to disk as they are made.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    capacity = len(recorded_steps(config.steps, frame_every))
    sink = EcosystemFrameSink.create(
        output_dir / "frames", config, capacity=capacity, overwrite=True
    )
    try:
        result = run_ecosystem(config, frame_every=frame_every, sink=sink)
    finally:
        index = sink.close()
    if result.manifest is None:
        raise RuntimeError("the ecosystem engine returned no manifest")
    store = FrameStore.open(output_dir / "frames")
    try:
        viewer = write_viewer(result, output_dir / "viewer.html", frames=StoredFrames(store))
    finally:
        store.close()  # release the memory maps, which Windows requires before files move
    manifest = result.manifest.write(output_dir / "manifest.json")
    publish_run(output_dir, result.manifest)
    return result, capacity, (index, viewer, manifest)


legacy_ecosystem_v1 = LegacyEcosystem()
