"""Recorded frames, streamed to disk as standard NumPy files.

A run can record hundreds of frames, each holding several fields over the whole
grid. Keeping them all in memory ties memory to run length, and writing them as
JSON multiplies their size several times over. A frame store instead gives
each field one preallocated ``.npy`` file with a leading frame axis, written
in place through a memory map. A run then holds one frame at a time, and any
tool that reads NumPy files can read the store without MARSE.

``index.json`` beside the arrays records the format, how many frames were
written, their times and steps, each field's shape and type, and metadata such
as the names that label each axis. The store knows nothing about any engine:
the ecosystem engine (:mod:`marse.ecosystem.framestore`) and the reaction-
transport engine (:mod:`marse.core.reactive_transport`) both write through it.

Frames are for looking at a run. The exact result is the final state recorded
in the manifest and reproduced by ``marse replay``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = ["FORMAT", "FORMAT_VERSION", "FrameStore", "FrameStoreError"]

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
