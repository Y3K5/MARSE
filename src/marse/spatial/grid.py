"""A box of cubic voxels over a flat substratum, in one, two or three dimensions.

The scene is the standard one for a biofilm: an impermeable surface (a tooth,
a pipe wall, a slide) at the bottom, colonies on it, and liquid above whose
bulk concentrations are held at the top face (docs/networks.md, "Running a
network in space"). The last axis is height above the substratum and the
others are lateral:

- 1-D: ``(nz,)``, a single column through the depth;
- 2-D: ``(nx, nz)``, a vertical slice;
- 3-D: ``(nx, ny, nz)``.

Lateral faces are periodic, so the box is one tile of a surface that repeats
sideways. A missing lateral axis stands for a depth of one voxel: a 1-D column
has the footprint of one voxel, and areal quantities (per m² of substratum)
mean the same thing in every dimension.

Fields on the grid have shape ``(components, *shape)``. Every operator in
:mod:`marse.spatial.transport` loops over the axes, so one code path serves
all three dimensions, and the lower-dimensional cases are exact special cases
of the higher ones.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse._numeric import require

__all__ = ["Grid"]


@dataclass(frozen=True, slots=True)
class Grid:
    """Voxels per axis, the last axis being height, and the voxel edge in micrometres."""

    shape: tuple[int, ...]
    voxel_um: float

    def __post_init__(self) -> None:
        require(1 <= len(self.shape) <= 3, "a grid has one, two or three axes")
        require(
            all(isinstance(n, int) and n >= 1 for n in self.shape),
            "voxel counts must be positive integers",
        )
        require(self.shape[-1] >= 2, "the height axis needs at least two voxels")
        require(math.isfinite(self.voxel_um) and self.voxel_um > 0, "voxel_um must be positive")

    @property
    def dimensions(self) -> int:
        return len(self.shape)

    @property
    def height_axis(self) -> int:
        return len(self.shape) - 1

    @property
    def lateral_axes(self) -> tuple[int, ...]:
        return tuple(range(len(self.shape) - 1))

    @property
    def voxels(self) -> int:
        return math.prod(self.shape)

    @property
    def size_um(self) -> tuple[float, ...]:
        return tuple(n * self.voxel_um for n in self.shape)

    @property
    def voxel_volume_um3(self) -> float:
        return self.voxel_um**3

    @property
    def footprint_um2(self) -> float:
        """Area of substratum the box covers; a missing lateral axis counts one voxel deep."""
        lateral = [self.shape[a] for a in self.lateral_axes]
        lateral += [1] * (2 - len(lateral))
        return math.prod(lateral) * self.voxel_um**2

    def centers_um(self, axis: int) -> NDArray[np.float64]:
        """Voxel centres along ``axis``; height is measured from the substratum."""
        return (np.arange(self.shape[axis]) + 0.5) * self.voxel_um

    def heights_um(self) -> NDArray[np.float64]:
        return self.centers_um(self.height_axis)
