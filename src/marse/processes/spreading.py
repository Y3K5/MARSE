"""Plaque that spreads as it grows, wears at its surface, and is brushed off, in a column.

Biomass does not diffuse; it is displaced (docs/theory.md, section 6.1). In a
column of voxels over the substratum, the components that occupy space each
have a packing concentration P, their ``density_mol_per_m3``: the
concentration at which they alone would fill a voxel. A domain spreads this
way when its ``spreading`` names the ``packed`` mechanism. A voxel's solid
fraction is

    phi = sum over the occupying components of X / P.

Growth raises phi above one. After every step of the integration the solid is
remapped so that it packs the column from the substratum up: the solid of old
voxel k lies between the cumulative volumes V[k - 1] and V[k], in voxels, and
new voxel j takes whatever lies between j and j + 1. Every component that
moves with the solid, the occupying ones and those it carries (such as the
fixed buffer of cell walls), goes with it, in the proportions of the voxel it
came from. That is the displacement of Wanner and Gujer (1986) on a fixed grid:
the velocity of the solid at a height is the growth of the solid below it.
The front is sharp, full voxels then at most one partly filled one, and where
growth slows below a voxel's capacity the column contracts. Dissolved
components stay where they are.

Three things take solid off the top:

- **The plaque's maximum height.** Whatever is pushed above it is detached:
  in the mouth, the film's flow carries it away.
- **Wear,** at a stated velocity: friction of the tongue and cheeks removes
  the surface. It is applied in two halves, before and after each step
  (Strang splitting), so it is second order in the step like the integration.
- **Removal** of a fraction of the plaque from its surface down, as brushing
  and flossing do.

The remap integrates each component's amount as a function of the solid volume
below it, a non-decreasing piecewise-linear function, and takes differences of
it between the new voxels' bounds. So every piece is non-negative, the pieces
and what is detached add up to what there was, to rounding, and a run replays
bit for bit. Spreading runs in a column only: how biomass spreads in two and
three dimensions changes conclusions, so that choice belongs to Stage 2d,
which compares mechanisms (docs/modeling-landscape.md, section 2).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

__all__ = ["SPREADING_VERSION", "Spreading", "remap"]

SPREADING_VERSION = "displacement_1d_v1"
"""The spreading mechanism, recorded in the manifest of every run that uses it."""

type Field = NDArray[np.float64]


def remap(
    column: Field,
    occupying: NDArray[np.intp],
    packing_mol_per_m3: NDArray[np.float64],
    moving: NDArray[np.intp],
    keep_voxels: float,
) -> tuple[Field, NDArray[np.float64]]:
    """Pack a column's solid from the substratum up, and keep at most ``keep_voxels`` of it.

    ``column`` has shape (components, voxels), in concentrations; ``moving``
    lists every component that moves with the solid. A voxel that holds no
    solid keeps what it holds. Returns the new column and what lies above the
    cut, per component, as concentration summed over voxels.
    """
    voxels = column.shape[1]
    phi = np.maximum((column[occupying] / packing_mol_per_m3[:, None]).sum(axis=0), 0.0)
    knots = np.concatenate(([0.0], np.cumsum(phi)))
    cut = min(float(keep_voxels), float(knots[-1]), float(voxels))
    # In a voxel without solid, nothing carries what the voxel holds: it stays.
    empty = phi == 0.0
    amounts = column[moving]
    carried = np.where(empty, 0.0, amounts)
    below = np.concatenate((np.zeros((amounts.shape[0], 1)), np.cumsum(carried, axis=1)), axis=1)
    bounds = np.minimum(np.arange(voxels + 1, dtype=float), max(cut, 0.0))
    at = np.empty((amounts.shape[0], voxels + 1))
    for i in range(amounts.shape[0]):
        at[i] = np.interp(bounds, knots, below[i], right=below[i, -1])
    pieces = np.maximum(np.diff(at, axis=1), 0.0)
    new = column.copy()
    new[moving] = pieces + np.where(empty, amounts, 0.0)
    detached = np.zeros(column.shape[0])
    detached[moving] = np.maximum(below[:, -1] - pieces.sum(axis=1), 0.0)
    return new, detached


@dataclass(frozen=True, slots=True)
class Spreading:
    """The solid phase of a column: what fills it, what it carries, and what takes it off.

    ``occupying`` and ``packing_mol_per_m3`` give the components that fill
    space and their packing concentrations; ``moving`` adds the components
    carried with them. The solid is at most ``maximum_voxels`` high, and wears
    at ``wear_um_per_h``.
    """

    occupying: NDArray[np.intp]
    packing_mol_per_m3: NDArray[np.float64]
    moving: NDArray[np.intp]
    maximum_voxels: float
    voxel_um: float
    wear_um_per_h: float = 0.0
    _wear_voxels_per_h: float = field(init=False, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "_wear_voxels_per_h", self.wear_um_per_h / self.voxel_um)

    def solid_fraction(self, box: Field) -> NDArray[np.float64]:
        """phi in every voxel, shaped like one component of ``box``."""
        column = box.reshape(box.shape[0], -1)
        phi = (column[self.occupying] / self.packing_mol_per_m3[:, None]).sum(axis=0)
        return phi.reshape(box.shape[1:])

    def height_um(self, box: Field) -> float:
        """The plaque's thickness: its total solid volume, as a height."""
        return float(self.solid_fraction(box).sum()) * self.voxel_um

    def project(self, box: Field, hours: float) -> tuple[Field, NDArray[np.float64]]:
        """Wear ``hours`` of the surface away, pack the rest, and detach what passes the top.

        Returns the new box and what left it, per component, as concentration
        summed over voxels.
        """
        column = box.reshape(box.shape[0], -1)
        height = float(self.solid_fraction(box).sum())
        keep = min(self.maximum_voxels, height - min(height, self._wear_voxels_per_h * hours))
        column, detached = remap(column, self.occupying, self.packing_mol_per_m3, self.moving, keep)
        return column.reshape(box.shape), detached

    def remove(self, box: Field, fraction: float) -> tuple[Field, NDArray[np.float64]]:
        """Take ``fraction`` of the plaque off from its surface down, as a brush does."""
        column = box.reshape(box.shape[0], -1)
        height = float(self.solid_fraction(box).sum())
        keep = min(self.maximum_voxels, height * (1.0 - fraction))
        column, removed = remap(column, self.occupying, self.packing_mol_per_m3, self.moving, keep)
        return column.reshape(box.shape), removed
