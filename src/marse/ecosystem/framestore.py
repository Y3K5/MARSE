"""Recorded frames, streamed to disk as standard NumPy files.

A run can record hundreds of frames, each holding several fields over the whole
grid. Keeping them all in memory ties memory to run length, and writing them as
JSON multiplies their size several times over. A frame store instead gives
each field one preallocated ``.npy`` file with a leading frame axis, written
in place through a memory map. A run then holds one frame at a time, and any
tool that reads NumPy files can read the store without MARSE.

``index.json`` beside the arrays records the format, how many frames were
written, their times and steps, each field's shape and type, and metadata such
as the names that label each axis. :class:`FrameStore` lives in
:mod:`marse.core.framestore` and knows nothing about ecosystems; this module
maps ecosystem frames onto it and re-exports it.

Frames are for looking at a run. The exact result is the final state recorded
in the manifest and reproduced by ``marse replay``. The ecosystem adapter
therefore stores fields in single precision by default.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from pathlib import Path

import numpy as np

from marse.core.framestore import FORMAT, FORMAT_VERSION, FrameStore, FrameStoreError
from marse.core.provenance import config_checksum
from marse.ecosystem.model import EcosystemConfig, EcosystemError, EcosystemFrame

__all__ = [
    "FORMAT",
    "FORMAT_VERSION",
    "EcosystemFrameSink",
    "FrameStore",
    "FrameStoreError",
    "StoredFrames",
]


# --- the ecosystem adapter -----------------------------------------------------


def _limiting_legend(config: EcosystemConfig) -> tuple[str, ...]:
    """Every name a limiting-factor map can hold, so the maps can be stored as codes."""
    names = {""}
    for species in config.species:
        for capability in species.capabilities:
            names.add(capability.substrate)
            if capability.temperature_c is not None:
                names.add("temperature")
            if capability.ph is not None:
                names.add("ph")
            if capability.oxygen_half_saturation is not None:
                names.add("oxygen")
    return tuple(sorted(names))


class EcosystemFrameSink:
    """Streams ecosystem frames into a :class:`FrameStore`, for ``run(..., sink=...)``."""

    def __init__(self, store: FrameStore, legend: tuple[str, ...]) -> None:
        self.store = store
        self._legend = np.asarray(legend)

    @classmethod
    def create(
        cls,
        directory: str | Path,
        config: EcosystemConfig,
        *,
        capacity: int,
        precision: str = "float32",
        overwrite: bool = False,
    ) -> EcosystemFrameSink:
        if precision not in ("float32", "float64"):
            raise EcosystemError(f"precision must be float32 or float64, got {precision!r}")
        grid = (config.height, config.width)
        s, n = len(config.species), len(config.nutrients)
        c, a = len(config.conditions), len(config.additives)
        real = f"<f{4 if precision == 'float32' else 8}"
        legend = _limiting_legend(config)
        fields = {
            "biomass": ((s, *grid), real),
            "nutrients": ((n, *grid), real),
            "conditions": ((c, *grid), real),
            "additives": ((a, *grid), real),
            "mutations": ((s, *grid), "<u1"),
            "phenotype_indices": ((s, *grid), "<u2"),
            "niche_rates": ((s, *grid), real),
            "limiting_factor_codes": ((s, *grid), "<u1"),
        }
        metadata = {
            "experiment_id": config.experiment_id,
            "config_sha256": config_checksum(config),
            "width": config.width,
            "height": config.height,
            "cell_size_um": config.cell_size_um,
            "carrying_capacity": config.carrying_capacity,
            "species": [x.name for x in config.species],
            "nutrients": [x.name for x in config.nutrients],
            "conditions": [x.name for x in config.conditions],
            "additives": [x.name for x in config.additives],
            "limiting_factor_legend": list(legend),
            "precision": precision,
        }
        store = FrameStore.create(
            directory, capacity=capacity, fields=fields, metadata=metadata, overwrite=overwrite
        )
        return cls(store, legend)

    def write(self, index: int, step: int, frame: EcosystemFrame) -> None:
        limiting = np.stack(frame.niche_limiting_factors)
        codes = np.searchsorted(self._legend, limiting)
        codes = np.minimum(codes, len(self._legend) - 1)
        if not np.array_equal(self._legend[codes], limiting):
            raise EcosystemError("a limiting factor outside the configured names was recorded")
        written = self.store.append(
            time_h=frame.time_h,
            step=step,
            arrays={
                "biomass": frame.biomass,
                "nutrients": frame.nutrients,
                "conditions": frame.conditions,
                "additives": frame.additives,
                "mutations": frame.mutations,
                "phenotype_indices": (
                    frame.phenotype_indices
                    if frame.phenotype_indices is not None
                    else np.zeros_like(frame.mutations)
                ),
                "niche_rates": np.stack(frame.niche_rates),
                "limiting_factor_codes": codes,
            },
        )
        if written != index:
            raise EcosystemError(f"frame {index} was stored at position {written}")

    def close(self) -> Path:
        return self.store.close()


class StoredFrames(Sequence[EcosystemFrame]):
    """Frames read back from a store on demand, one at a time.

    Behaves like ``result.frames`` for the viewer and for analysis, without
    loading the whole run into memory.
    """

    def __init__(self, store: FrameStore) -> None:
        self.store = store
        self._legend = np.asarray(store.metadata["limiting_factor_legend"])

    def __len__(self) -> int:
        return len(self.store)

    def __getitem__(self, index):  # type: ignore[override]
        if isinstance(index, slice):
            return [self[i] for i in range(*index.indices(len(self)))]
        if index < 0:
            index += len(self)
        fields = self.store.read(index)
        return EcosystemFrame(
            self.store.times_h[index],
            fields["biomass"].astype(float),
            fields["nutrients"].astype(float),
            fields["conditions"].astype(float),
            fields["additives"].astype(float),
            fields["mutations"].astype(np.int64),
            tuple(fields["niche_rates"].astype(float)),
            tuple(self._legend[fields["limiting_factor_codes"]]),
            fields["phenotype_indices"].astype(np.int64),
        )

    def __iter__(self) -> Iterator[EcosystemFrame]:
        for index in range(len(self)):
            yield self[index]
