"""The salivary film over a site of plaque: how fast fresh saliva replaces it.

Saliva covers the teeth as a film about 0.1 mm thick (Collins and Dawes 1987)
that moves over them at 0.8 to 7.6 mm per minute, depending on the site
(Dawes et al. 1989). Clearance from plaque into the film takes longer the more
plaque the film has already crossed, and less long the faster it moves (Dawes
1989). A film of thickness delta and mean velocity u_bar, with its surface
open to the air, moves at

    u(zeta) = 1.5 u_bar (2 zeta - zeta^2),   zeta = height in the film / delta,

zero at the plaque and fastest at the surface; its shear at the plaque is
3 u_bar / delta, as the dental scene of increment E1 uses (docs/theory.md,
section 6.4). Over a site that the film reaches after crossing a length l of
plaque, each layer of the film is replaced at the rate u(zeta) / l
(docs/theory.md, section 4.8): a well-mixed film renewed in plug flow, with
the same steady state as the film leaving the plaque.

Over a plaque that spreads (:mod:`marse.biofilm.spreading`), the film rides on
the plaque's surface, wherever it is: zeta is measured from the surface, and
each voxel's liquid, the part of it above the surface, is renewed. Liquid
more than a film's thickness above the surface, which the box leaves when
plaque is brushed off below its maximum height, moves with the film's
surface (docs/theory.md, section 4.8).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from marse.schemas.domain import Film
from marse.spatial.grid import Grid

__all__ = ["film_layers", "film_share", "liquid_share", "renewal_over", "renewal_per_h"]

_MM_PER_MIN_TO_UM_PER_H = 1000.0 * 60.0


def film_layers(grid: Grid, film: Film) -> int:
    """How many voxels deep the film is: the top layers of the box."""
    return round(film.thickness_um / grid.voxel_um)


def renewal_per_h(grid: Grid, film: Film) -> NDArray[np.float64]:
    """The rate at which fresh saliva replaces each voxel of the box, per hour.

    Zero below the film; in the film's voxels, u(zeta) / l at their centres.
    """
    layers = film_layers(grid, film)
    zeta = (np.arange(layers) + 0.5) / layers
    speed = 1.5 * film.velocity_mm_per_min * _MM_PER_MIN_TO_UM_PER_H * (2 * zeta - zeta**2)
    rate = np.zeros(grid.shape)
    rate[..., -layers:] = speed / (film.plaque_length_mm * 1000.0)
    return rate


def liquid_share(grid: Grid, surface_um: float) -> NDArray[np.float64]:
    """The share of each voxel above the plaque's surface, ``surface_um`` up: its liquid."""
    h = grid.voxel_um
    tops = (np.arange(grid.shape[-1]) + 1.0) * h
    share = np.clip((tops - surface_um) / h, 0.0, 1.0)
    return np.broadcast_to(share, grid.shape).copy()


def renewal_over(grid: Grid, film: Film, surface_um: float) -> NDArray[np.float64]:
    """The renewal of each voxel, per hour, under a film on a plaque ``surface_um`` thick.

    Each voxel's liquid is renewed at u(zeta) / l at the liquid's centre, with
    zeta its height above the surface over the film's thickness, at most one;
    the voxel's renewal is that times its share of liquid. On a plaque as high
    as the film's underside, this is :func:`renewal_per_h`, to rounding.
    """
    h = grid.voxel_um
    share = np.clip((np.arange(grid.shape[-1]) + 1.0 - surface_um / h), 0.0, 1.0)
    centres = (np.arange(grid.shape[-1]) + 1.0 - 0.5 * share) * h
    zeta = np.clip((centres - surface_um) / film.thickness_um, 0.0, 1.0)
    speed = 1.5 * film.velocity_mm_per_min * _MM_PER_MIN_TO_UM_PER_H * (2 * zeta - zeta**2)
    rate = share * speed / (film.plaque_length_mm * 1000.0)
    return np.broadcast_to(rate, grid.shape).copy()


def film_share(grid: Grid, film: Film, surface_um: float) -> NDArray[np.float64]:
    """The share of each voxel in the film's thickness above the plaque's surface.

    Over the column it adds up to the film's thickness, in voxels, as long as
    the surface is no higher than the film's underside.
    """
    h = grid.voxel_um
    bottoms = np.arange(grid.shape[-1]) * h
    overlap = np.minimum(bottoms + h, surface_um + film.thickness_um) - np.maximum(
        bottoms, surface_um
    )
    return np.broadcast_to(np.clip(overlap / h, 0.0, 1.0), grid.shape).copy()
