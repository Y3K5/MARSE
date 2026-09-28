"""Finite-volume diffusion on a grid of voxels, conservative by construction.

Concentrations ``c`` have shape ``(components, *shape)`` on a
:class:`~marse.spatial.grid.Grid`. Along every axis, the flux through the face
between two voxels is computed once, from their difference,

    F = -D (c_upper - c_lower) / h,

and the same value leaves one voxel and enters the other. Whatever a voxel
loses, a neighbour gains, so the operator can neither create nor destroy
matter; only the top face exchanges with the bulk liquid, and that exchange
is what the ledger counts as imports and exports (docs/theory.md, sections
4.7 and 9.5).

Boundaries (docs/networks.md, "Running a network in space"):

- **lateral axes are periodic**: the box is one tile of a repeating surface;
- **the substratum** (below the first voxel in height) is impermeable;
- **the top face** is held at the bulk-liquid concentrations. The voxel centre
  lies half a voxel below it, so the flux there is ``-D (bulk - c) / (h / 2)``.
  Under a salivary film, whose surface is open to the air, the top face is
  closed instead (``closed_top``), and the film exchanges with the mouth
  through its voxels (:mod:`marse.core.reservoir`).

Components whose diffusivity is zero (biomass) do not move.

Spacing is per axis. The model grid is cubic, but the coarse levels of the
multigrid solver (:mod:`marse.spatial.multigrid`) need not be, so every
function here takes the spacing of each axis.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from marse._numeric import require
from marse.spatial.grid import Grid

__all__ = ["Diffusion", "diagonal", "divergence", "face_fluxes"]


def _per_component(values: NDArray[np.float64], dims: int) -> NDArray[np.float64]:
    return values.reshape((-1,) + (1,) * dims)


def face_fluxes(
    c: NDArray[np.float64],
    diffusivity: NDArray[np.float64],
    spacing: tuple[float, ...],
    bulk: NDArray[np.float64] | None,
    *,
    closed_top: bool = False,
) -> tuple[NDArray[np.float64], ...]:
    """Flux through the upper face of every voxel along each axis.

    Returns one array per axis, each shaped like ``c``, positive in the
    direction of increasing index. On the height axis the last entry is the
    flux through the top face into the bulk liquid; ``bulk=None`` holds the
    bulk at zero, which is what the correction equations of Newton's method
    need, and ``closed_top`` makes it zero. The substratum face carries
    nothing and is not stored.
    """
    dims = c.ndim - 1
    d = _per_component(diffusivity, dims)
    fluxes = []
    for axis in range(dims - 1):  # lateral, periodic
        fluxes.append(-d * (np.roll(c, -1, axis=axis + 1) - c) / spacing[axis])
    h = spacing[-1]
    height = np.empty_like(c)
    height[..., :-1] = -d * (c[..., 1:] - c[..., :-1]) / h
    if closed_top:
        height[..., -1] = 0.0
    else:
        top = c[..., -1]
        outside = 0.0 if bulk is None else _per_component(bulk, dims - 1)
        height[..., -1] = -_per_component(diffusivity, dims - 1) * (outside - top) / (0.5 * h)
    fluxes.append(height)
    return tuple(fluxes)


def divergence(
    fluxes: tuple[NDArray[np.float64], ...], spacing: tuple[float, ...]
) -> NDArray[np.float64]:
    """The rate of change each voxel sees: what enters through its lower faces less what leaves."""
    dims = len(fluxes)
    rate = np.zeros_like(fluxes[0])
    for axis in range(dims - 1):
        upper = fluxes[axis]
        rate += (np.roll(upper, 1, axis=axis + 1) - upper) / spacing[axis]
    upper = fluxes[-1]
    rate[..., 1:] += upper[..., :-1] / spacing[-1]
    rate -= upper / spacing[-1]
    return rate


def diagonal(
    diffusivity: NDArray[np.float64],
    spacing: tuple[float, ...],
    shape: tuple[int, ...],
    *,
    closed_top: bool = False,
) -> NDArray[np.float64]:
    """How each voxel's own rate depends on its own concentration: the diagonal of the operator."""
    dims = len(shape)
    d = _per_component(diffusivity, dims)
    diag = np.zeros((diffusivity.size, *shape))
    for axis in range(dims - 1):
        if shape[axis] >= 2:  # a periodic axis of one voxel exchanges with itself: nothing
            diag -= 2.0 * d / spacing[axis] ** 2
    h2 = spacing[-1] ** 2
    faces_below = np.ones(shape[-1])
    faces_below[0] = 0.0  # the substratum
    faces_above = np.ones(shape[-1])
    faces_above[-1] = 0.0 if closed_top else 2.0  # the top face, half a voxel away
    diag -= d * (faces_below + faces_above) / h2
    return diag


@dataclass(frozen=True, slots=True)
class Diffusion:
    """Diffusion of several components on a grid, with the bulk liquid held above it.

    ``diffusivity_um2_per_h`` and ``bulk_mol_per_m3`` hold one value per
    component. Amounts are concentration times volume, mol/m3 x um3, which is
    attomoles (1e-18 mol). With ``closed_top`` nothing crosses the top face,
    and the bulk is not used.
    """

    grid: Grid
    diffusivity_um2_per_h: NDArray[np.float64]
    bulk_mol_per_m3: NDArray[np.float64]
    closed_top: bool = False
    spacing: tuple[float, ...] = field(init=False)

    def __post_init__(self) -> None:
        d = np.asarray(self.diffusivity_um2_per_h, dtype=float)
        bulk = np.asarray(self.bulk_mol_per_m3, dtype=float)
        require(
            d.ndim == 1 and bulk.shape == d.shape, "one diffusivity and bulk value per component"
        )
        require(
            bool(np.all(np.isfinite(d)) and np.all(d >= 0)), "diffusivities must not be negative"
        )
        require(
            bool(np.all(np.isfinite(bulk)) and np.all(bulk >= 0)),
            "bulk concentrations must not be negative",
        )
        object.__setattr__(self, "diffusivity_um2_per_h", d)
        object.__setattr__(self, "bulk_mol_per_m3", bulk)
        object.__setattr__(self, "spacing", (self.grid.voxel_um,) * self.grid.dimensions)

    @property
    def components(self) -> int:
        return self.diffusivity_um2_per_h.size

    def fluxes(
        self, c: NDArray[np.float64], *, homogeneous: bool = False
    ) -> tuple[NDArray[np.float64], ...]:
        return face_fluxes(
            c,
            self.diffusivity_um2_per_h,
            self.spacing,
            None if homogeneous else self.bulk_mol_per_m3,
            closed_top=self.closed_top,
        )

    def rate(self, c: NDArray[np.float64], *, homogeneous: bool = False) -> NDArray[np.float64]:
        """dc/dt from diffusion alone, mol/m3 per hour."""
        return divergence(self.fluxes(c, homogeneous=homogeneous), self.spacing)

    def diagonal(self) -> NDArray[np.float64]:
        return diagonal(
            self.diffusivity_um2_per_h, self.spacing, self.grid.shape, closed_top=self.closed_top
        )

    def top_flux(self, c: NDArray[np.float64]) -> NDArray[np.float64]:
        """Flux into the box through each top face, shape (components, *lateral), mol/m3 x um/h."""
        return -self.fluxes(c)[-1][..., -1]

    def import_rate(self, c: NDArray[np.float64]) -> NDArray[np.float64]:
        """Amount entering through the top per hour, per component, in attomoles per hour."""
        top = self.top_flux(c)
        face_area = self.grid.voxel_um**2
        return top.reshape(top.shape[0], -1).sum(axis=1) * face_area

    def amounts(self, c: NDArray[np.float64]) -> NDArray[np.float64]:
        """Amount of each component in the box, attomoles."""
        return c.reshape(c.shape[0], -1).sum(axis=1) * self.grid.voxel_volume_um3

    def explicit_step_limit_h(self) -> float:
        """The longest stable explicit step, h^2 / (2 d D_max), which implicit steps avoid."""
        dmax = float(self.diffusivity_um2_per_h.max())
        if dmax == 0.0:
            return math.inf
        return self.grid.voxel_um**2 / (2 * self.grid.dimensions * dmax)
