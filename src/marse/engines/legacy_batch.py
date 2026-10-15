"""The version-1 well-mixed batch kernel.

This is the engine ``marse run`` uses for an experiment that is not schema
version 2 and has no biofilm block.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from marse.core.config import ConfigError, ExperimentConfig, experiment_from_dict, load_experiment
from marse.core.provenance import Manifest
from marse.core.simulation import run as simulate
from marse.engines.base import Completed
from marse.io.results import write_run_files


class LegacyBatch:
    """Well-mixed batch runs of the version-1 experiment schema."""

    name = "legacy_batch"
    kinds: tuple[str, ...] = ("batch",)
    unverified = False
    writes_during_replay = False

    def load(self, path: Path) -> ExperimentConfig:
        config = load_experiment(path)
        if config.biofilm is not None:
            raise ConfigError(
                "engine 'legacy_batch' runs a well-mixed batch; this document sets biofilm. "
                "Use engine 'steady_biofilm'."
            )
        return config

    def validate(self, document: dict[str, Any]) -> ExperimentConfig:
        config = experiment_from_dict(document)
        if config.biofilm is not None:
            raise ConfigError(
                "engine 'legacy_batch' runs a well-mixed batch; this document sets biofilm. "
                "Use engine 'steady_biofilm'."
            )
        return config

    def default_output_dir(self, plan: Any, source: Path) -> Path:
        return source.parent / "runs" / Manifest.run_id_for(plan)

    def run(
        self,
        plan: Any,
        output_dir: Path,
        *,
        frame_interval_h: float | None = None,
    ) -> Completed:
        if not isinstance(plan, ExperimentConfig):
            raise TypeError("legacy_batch runs an ExperimentConfig")
        result = simulate(plan)
        return Completed(result, write_run_files(result, output_dir))

    def replay(self, manifest: Manifest, output_dir: Path | None) -> Completed:
        config = _batch(manifest)
        return Completed(simulate(config))

    def write_outputs(self, result: Any, output_dir: Path) -> tuple[Path, ...]:
        return write_run_files(result, output_dir)


def _batch(manifest: Manifest) -> ExperimentConfig:
    config = manifest.experiment()
    if not isinstance(config, ExperimentConfig) or config.kind != "batch":
        raise ValueError(
            f"manifest declares {manifest.kind!r} but its configuration is not a batch experiment"
        )
    return config


legacy_batch = LegacyBatch()
