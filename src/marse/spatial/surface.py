"""The substratum's materials: which surface every face of the bottom layer is.

A scene may pattern its substratum, putting enamel, titanium, zirconia and
acrylic side by side under one liquid, or bare and saliva-coated glass in one
flow chamber (docs/environments.md). Each face of the substratum, one below
every column of voxels, is made of exactly one material. Rectangular patches,
stated in micrometres, assign faces by where their centres lie, so a face is
never split between two materials and none is left without one.

The face is the unit a material, and with it adhesion, attaches to. Stage E3
gives the solid grains of a soil faces of their own, in every direction; the
substratum is the first case of the same idea.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse.spatial.grid import Grid

__all__ = ["Patch", "Substratum"]


@dataclass(frozen=True, slots=True)
class Patch:
    """A rectangle of the substratum made of one material.

    ``region_um`` holds a lower and an upper bound for each lateral axis in
    turn: ``(x0, x1, y0, y1)`` in 3-D, ``(x0, x1)`` in 2-D, and nothing in 1-D,
    where a column has a single face.
    """

    material: str
    region_um: tuple[float, ...]


class Substratum:
    """The material of every face of the substratum, from non-overlapping patches."""

    def __init__(self, grid: Grid, patches: Sequence[Patch]) -> None:
        if not patches:
            raise ValueError("the substratum needs at least one patch")
        lateral = grid.lateral_axes
        shape = tuple(grid.shape[a] for a in lateral)
        owner = np.full(shape, -1, dtype=np.intp)
        materials: list[str] = []
        for number, patch in enumerate(patches):
            if len(patch.region_um) != 2 * len(lateral):
                raise ValueError(
                    f"patch {number} ({patch.material}): a {grid.dimensions}-D box needs "
                    f"{2 * len(lateral)} bounds (a lower and an upper one per lateral axis), "
                    f"got {len(patch.region_um)}"
                )
            inside = np.ones(shape, dtype=bool)
            for k, axis in enumerate(lateral):
                low, high = patch.region_um[2 * k], patch.region_um[2 * k + 1]
                size = grid.size_um[axis]
                if not 0 <= low < high <= size:
                    raise ValueError(
                        f"patch {number} ({patch.material}): bounds {low:g} to {high:g} um must "
                        f"satisfy 0 <= lower < upper <= {size:g}"
                    )
                centres = grid.centers_um(axis)
                along = (centres >= low) & (centres < high)
                inside &= along.reshape([-1 if j == k else 1 for j in range(len(lateral))])
            if not inside.any():
                raise ValueError(
                    f"patch {number} ({patch.material}): covers no face; faces are assigned by "
                    "their centres"
                )
            taken = owner[inside]
            if (taken >= 0).any():
                other = patches[int(taken[taken >= 0][0])].material
                raise ValueError(
                    f"patch {number} ({patch.material}) overlaps an earlier patch ({other})"
                )
            owner[inside] = number
            if patch.material not in materials:
                materials.append(patch.material)
        uncovered = int((owner < 0).sum())
        if uncovered:
            raise ValueError(
                f"{uncovered} face(s) of the substratum belong to no patch; the patches must "
                "cover it"
            )
        material_of_patch = np.array([materials.index(p.material) for p in patches], dtype=np.intp)
        self.grid = grid
        self.patches = tuple(patches)
        self.materials = tuple(materials)
        self.index: NDArray[np.intp] = material_of_patch[owner]

    @property
    def face_area_um2(self) -> float:
        """Each face's area; a missing lateral axis counts one voxel deep, as in the footprint."""
        return self.grid.voxel_um**2

    def mask(self, material: str) -> NDArray[np.bool_]:
        """Which faces are made of ``material``."""
        return self.index == self.materials.index(material)

    def area_um2(self, material: str) -> float:
        return float(self.mask(material).sum()) * self.face_area_um2

    def per_face(self, values: Mapping[str, float]) -> NDArray[np.float64]:
        """A value for each face, from a value for each material."""
        missing = [m for m in self.materials if m not in values]
        if missing:
            raise ValueError(f"no value for material(s) {', '.join(missing)}")
        table = np.array([float(values[m]) for m in self.materials])
        return table[self.index]
