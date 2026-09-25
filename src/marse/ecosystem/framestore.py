"""Recorded frames, streamed to disk as standard NumPy files.

A run can record hundreds of frames, each holding several fields over the whole
grid. Keeping them all in memory ties memory to run length, and writing them as
JSON multiplies their size several times over. A frame store instead gives
each field one preallocated ``.npy`` file with a leading frame axis, written
in place through a memory map. A run then holds one frame at a time, and any
tool that reads NumPy files can read the store without MARSE.

``index.json`` beside the arrays records the format, how many frames were
written, their times and steps, each field's shape and type, and metadata such
as the names that label each axis. :class:`FrameStore` knows nothing about
ecosystems. The adapter at the end of this module maps ecosystem frames onto
it, so a later engine can reuse the store with different fields.

Frames are for looking at a run. The exact result is the final state recorded
in the manifest and reproduced by ``marse replay``. The ecosystem adapter
therefore stores fields in single precision by default.
"""

from __future__ import annotations

import json
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from marse.core.provenance import config_checksum
from marse.ecosystem.model import EcosystemConfig, EcosystemError, EcosystemFrame

__all__ = [
    "EcosystemFrameSink",
    "FrameStore",
    "FrameStoreError",
    "StoredFrames",
]

FORMAT = "marse-frames"
FORMAT_VERSION = 1


class FrameStoreError(ValueError):
    """A frame store is malformed, full, or used in the wrong mode."""


class FrameStore:
    """Named arrays with a leading frame axis, one ``.npy`` file per field.

    Create one with :meth:`create` and fill it with :meth:`append`, then call
    :meth:`close`, which writes ``index.json``. Read one with :meth:`open`.
    Fields with no elements (an ecosystem without additives, say) are recorded
    in the index but get no file, because an empty file cannot be mapped.
    """

    def __init__(
        self,
        directory: Path,
        fields: dict[str, tuple[tuple[int, ...], np.dtype]],
        arrays: dict[str, np.ndarray],
        capacity: int,
        metadata: dict[str, Any],
        times_h: list[float],
        steps: list[int],
        *,
        writable: bool,
    ) -> None:
        self.directory = directory
        self.fields = fields
        self.capacity = capacity
        self.metadata = metadata
        self._arrays = arrays
        self._times_h = times_h
        self._steps = steps
        self._writable = writable

    @classmethod
    def create(
        cls,
        directory: str | Path,
        *,
        capacity: int,
        fields: dict[str, tuple[tuple[int, ...], str]],
        metadata: dict[str, Any] | None = None,
        overwrite: bool = False,
    ) -> FrameStore:
        """Create an empty store. An existing store is refused unless ``overwrite``.

        Overwriting removes only the files the old store's index names, never
        anything else in the directory.
        """
        if capacity < 1:
            raise FrameStoreError("a frame store needs room for at least one frame")
        root = Path(directory)
        root.mkdir(parents=True, exist_ok=True)
        index_path = root / "index.json"
        if index_path.exists():
            if not overwrite:
                raise FrameStoreError(f"{root} already holds a frame store")
            old = json.loads(index_path.read_text(encoding="utf-8"))
            for name in old.get("fields", {}):
                (root / f"{name}.npy").unlink(missing_ok=True)
            index_path.unlink()
        layout = {name: (tuple(shape), np.dtype(dtype)) for name, (shape, dtype) in fields.items()}
        arrays = {
            name: np.lib.format.open_memmap(
                root / f"{name}.npy", mode="w+", dtype=dtype, shape=(capacity, *shape)
            )
            for name, (shape, dtype) in layout.items()
            if int(np.prod(shape)) > 0
        }
        return cls(root, layout, arrays, capacity, dict(metadata or {}), [], [], writable=True)

    def append(self, *, time_h: float, step: int, arrays: dict[str, ArrayLike]) -> int:
        """Write the next frame and return its index."""
        if not self._writable:
            raise FrameStoreError("this frame store was opened for reading")
        index = len(self._steps)
        if index >= self.capacity:
            raise FrameStoreError(f"the frame store is full ({self.capacity} frames)")
        missing = set(self.fields) - set(arrays)
        if missing:
            raise FrameStoreError(f"frame {index} is missing field(s) {sorted(missing)}")
        for name, (shape, _) in self.fields.items():
            values = np.asarray(arrays[name])
            if values.shape != shape:
                raise FrameStoreError(f"field {name!r}: expected shape {shape}, got {values.shape}")
            if name in self._arrays:
                self._arrays[name][index] = values
        self._times_h.append(float(time_h))
        self._steps.append(int(step))
        return index

    def close(self) -> Path:
        """Flush the arrays and write ``index.json``. Returns the index path."""
        if self._writable:
            for array in self._arrays.values():
                array.flush()
            index = {
                "format": FORMAT,
                "format_version": FORMAT_VERSION,
                "frame_count": len(self._steps),
                "capacity": self.capacity,
                # False when a run stopped early: the frames written are still valid.
                "complete": len(self._steps) == self.capacity,
                "times_h": self._times_h,
                "steps": self._steps,
                "fields": {
                    name: {"shape": list(shape), "dtype": dtype.str}
                    for name, (shape, dtype) in self.fields.items()
                },
                "metadata": self.metadata,
            }
            path = self.directory / "index.json"
            path.write_text(json.dumps(index, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            self._writable = False
        # Drop the maps so their files are released, which Windows requires.
        self._arrays = {}
        return self.directory / "index.json"

    @classmethod
    def open(cls, directory: str | Path) -> FrameStore:
        root = Path(directory)
        try:
            index = json.loads((root / "index.json").read_text(encoding="utf-8"))
        except FileNotFoundError as error:
            raise FrameStoreError(f"{root} has no index.json; was the store closed?") from error
        if index.get("format") != FORMAT or index.get("format_version") != FORMAT_VERSION:
            raise FrameStoreError(
                f"{root}: not a {FORMAT} store of version {FORMAT_VERSION} "
                f"(found {index.get('format')!r} version {index.get('format_version')!r})"
            )
        layout = {
            name: (tuple(spec["shape"]), np.dtype(spec["dtype"]))
            for name, spec in index["fields"].items()
        }
        arrays = {
            name: np.load(root / f"{name}.npy", mmap_mode="r")
            for name, (shape, _) in layout.items()
            if int(np.prod(shape)) > 0
        }
        count = int(index["frame_count"])
        return cls(
            root,
            layout,
            arrays,
            int(index["capacity"]),
            index["metadata"],
            [float(t) for t in index["times_h"][:count]],
            [int(s) for s in index["steps"][:count]],
            writable=False,
        )

    def __len__(self) -> int:
        return len(self._steps)

    @property
    def times_h(self) -> tuple[float, ...]:
        return tuple(self._times_h)

    @property
    def steps(self) -> tuple[int, ...]:
        return tuple(self._steps)

    def read(self, index: int) -> dict[str, NDArray]:
        """Every field of one frame, copied out of the maps."""
        if not 0 <= index < len(self):
            raise IndexError(f"frame {index} is outside 0..{len(self) - 1}")
        return {
            name: (
                np.array(self._arrays[name][index])
                if name in self._arrays
                else np.zeros(shape, dtype=dtype)
            )
            for name, (shape, dtype) in self.fields.items()
        }


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
