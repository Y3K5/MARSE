"""The one-dimensional steady biofilm solver.

A version-1 experiment with a ``biofilm`` block. The oxygen profile is checked
against analytical limits (validation case V3). The solver itself still lives
in :mod:`marse.core.simulation` and :mod:`marse.spatial.diffusion`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from marse.core.config import ConfigError, ExperimentConfig, experiment_from_dict, load_experiment
from marse.core.provenance import Manifest
from marse.core.simulation import run as simulate
from marse.engines.base import Completed
from marse.io.results import write_run_files


class SteadyBiofilm:
    """Steady reaction-diffusion through a biofilm of fixed biomass."""

    name = "steady_biofilm"
    kinds: tuple[str, ...] = ("biofilm_profile",)
    unverified = False
    writes_during_replay = False

    def load(self, path: Path) -> ExperimentConfig:
        config = load_experiment(path)
        if config.biofilm is None:
            raise ConfigError(
                "engine 'steady_biofilm' needs a biofilm block; this document has none. "
                "Use engine 'legacy_batch'."
            )
        return config

    def validate(self, document: dict[str, Any]) -> ExperimentConfig:
        config = experiment_from_dict(document)
        if config.biofilm is None:
            raise ConfigError(
                "engine 'steady_biofilm' needs a biofilm block; this document has none. "
                "Use engine 'legacy_batch'."
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
            raise TypeError("steady_biofilm runs an ExperimentConfig")
        result = simulate(plan)
        return Completed(result, write_run_files(result, output_dir))

    def replay(self, manifest: Manifest, output_dir: Path | None) -> Completed:
        config = manifest.experiment()
        if not isinstance(config, ExperimentConfig) or config.kind != "biofilm_profile":
            raise ValueError(
                f"manifest declares {manifest.kind!r} but its configuration "
                "is not a biofilm profile"
            )
        return Completed(simulate(config))

    def write_outputs(self, result: Any, output_dir: Path) -> tuple[Path, ...]:
        return write_run_files(result, output_dir)


steady_biofilm = SteadyBiofilm()
