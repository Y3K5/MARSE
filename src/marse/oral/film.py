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
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from marse.schemas.domain import Film
from marse.spatial.grid import Grid

__all__ = ["film_layers", "renewal_per_h"]

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
