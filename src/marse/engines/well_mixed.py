"""The version-2 reaction network in a closed, well-mixed box."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from marse.core.config import ConfigError
from marse.core.provenance import Manifest
from marse.core.well_mixed import run as run_well_mixed
from marse.engines.base import Completed
from marse.io.results import write_run_files
from marse.schemas.experiment import WellMixedConfig
from marse.schemas.experiment import experiment_from_dict as network_from_dict
from marse.schemas.experiment import load_experiment as load_network


class WellMixed:
    """Integrate a balanced reaction network with no spatial domain."""

    name = "well_mixed"
    kinds: tuple[str, ...] = ("well_mixed",)
    unverified = False
    writes_during_replay = False

    def load(self, path: Path) -> WellMixedConfig:
        config = load_network(path)
        if not isinstance(config, WellMixedConfig):
            raise ConfigError(
                "engine 'well_mixed' cannot run a document that sets domain; "
                "use engine 'reactive_transport'"
            )
        return config

    def validate(self, document: dict[str, Any]) -> WellMixedConfig:
        config = network_from_dict(document)
        if not isinstance(config, WellMixedConfig):
            raise ConfigError(
                "engine 'well_mixed' cannot run a document that sets domain; "
                "use engine 'reactive_transport'"
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
        if not isinstance(plan, WellMixedConfig):
            raise TypeError("well_mixed runs a WellMixedConfig")
        result = run_well_mixed(plan)
        return Completed(result, write_run_files(result, output_dir))

    def replay(self, manifest: Manifest, output_dir: Path | None) -> Completed:
        config = manifest.experiment()
        if not isinstance(config, WellMixedConfig):
            raise ValueError(
                f"manifest declares {manifest.kind!r} but its configuration is not well mixed"
            )
        return Completed(run_well_mixed(config))

    def write_outputs(self, result: Any, output_dir: Path) -> tuple[Path, ...]:
        return write_run_files(result, output_dir)


well_mixed = WellMixed()
