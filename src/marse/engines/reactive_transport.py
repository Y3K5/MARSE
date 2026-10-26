"""The version-2 reaction network in a box of voxels.

Running writes a frame store, VTK files and the totals the manifest records.
The numerical engine remains :mod:`marse.core.reactive_transport`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from marse.core.config import ConfigError
from marse.core.provenance import Manifest
from marse.core.reactive_transport import ReactiveTransportResult
from marse.core.reactive_transport import run as run_reactive_transport
from marse.engines.base import Completed
from marse.io.frames import FrameStore
from marse.io.results import write_run_files
from marse.model.experiment import ReactiveTransportConfig
from marse.model.experiment import experiment_from_dict as network_from_dict
from marse.model.experiment import load_experiment as load_network
from marse.processes.acid_base import ChargeBalance
from marse.spatial.vtk import write_pvd, write_vti


class ReactiveTransport:
    """Implicit reaction-transport on a one-, two- or three-dimensional grid."""

    name = "reactive_transport"
    kinds: tuple[str, ...] = ("reactive_transport",)
    unverified = False
    writes_during_replay = True

    def load(self, path: Path) -> ReactiveTransportConfig:
        config = load_network(path)
        if not isinstance(config, ReactiveTransportConfig):
            raise ConfigError(
                "engine 'reactive_transport' needs a domain; this document has none. "
                "Use engine 'well_mixed'."
            )
        return config

    def validate(self, document: dict[str, Any]) -> ReactiveTransportConfig:
        config = network_from_dict(document)
        if not isinstance(config, ReactiveTransportConfig):
            raise ConfigError(
                "engine 'reactive_transport' needs a domain; this document has none. "
                "Use engine 'well_mixed'."
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
        if not isinstance(plan, ReactiveTransportConfig):
            raise TypeError("reactive_transport runs a ReactiveTransportConfig")
        result, written = _run_in_space(plan, output_dir)
        return Completed(result, written)

    def replay(self, manifest: Manifest, output_dir: Path | None) -> Completed:
        config = manifest.experiment()
        if not isinstance(config, ReactiveTransportConfig):
            raise ValueError(
                f"manifest declares {manifest.kind!r} but its configuration is not in space"
            )
        if output_dir is None:
            return Completed(run_reactive_transport(config))
        return self.run(config, output_dir)

    def write_outputs(self, result: Any, output_dir: Path) -> tuple[Path, ...]:
        return write_run_files(result, output_dir)


def _run_in_space(
    config: ReactiveTransportConfig, output_dir: Path
) -> tuple[ReactiveTransportResult, tuple[Path, ...]]:
    """Run in space, streaming every recorded frame to a frame store and to VTK files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    grid = config.domain.grid
    names = config.network.component_names
    records = 1 + sum(
        1
        for step in range(1, config.steps + 1)
        if step % config.record_every == 0 or step == config.steps
    )
    balance = ChargeBalance.of(config.network) if config.network.has_ph else None
    layout = {f"{n}_mol_per_m3": (grid.shape, "<f4") for n in names}
    if balance is not None:
        layout["ph"] = (grid.shape, "<f4")
    store = FrameStore.create(
        output_dir / "frames",
        capacity=records,
        fields=layout,
        metadata={"voxel_um": grid.voxel_um, "axes": ["x", "y", "z"][-grid.dimensions :]},
        overwrite=True,
    )
    vtk = output_dir / "vtk"
    vtk.mkdir(exist_ok=True)
    series: list[tuple[float, str]] = []

    def record(index: int, time_h: float, state: np.ndarray) -> None:
        fields = {f"{n}_mol_per_m3": state[j] for j, n in enumerate(names)}
        if balance is not None:
            fields["ph"] = balance.ph(state)
        store.append(time_h=time_h, step=index, arrays=fields)
        name = f"frame_{index:04d}.vti"
        write_vti(vtk / name, grid, fields)
        series.append((time_h, name))

    try:
        result = run_reactive_transport(config, frames=record)
    finally:
        index = store.close()
    pvd = write_pvd(vtk / "run.pvd", series)
    totals, manifest = write_run_files(result, output_dir)
    written: list[Path] = [totals]
    if result.surface is not None:
        written.append(result.write_surface(output_dir / "surface.csv"))
    if result.ph is not None:
        written.append(result.write_ph(output_dir / "ph.csv"))
    if result.mouth is not None:
        written.append(result.write_mouth(output_dir / "mouth.csv"))
    if result.plaque is not None:
        written.append(result.write_plaque(output_dir / "plaque.csv"))
    if result.structure is not None:
        written.append(result.write_structure(output_dir / "structure.csv"))
    return result, (*written, index, pvd, manifest)


reactive_transport = ReactiveTransport()
