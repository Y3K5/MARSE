"""Cells binding to the substratum: deposition, blocking, detachment and locking.

The spatial engine exchanges cells between the suspension in the bulk liquid
and the faces of the substratum (docs/theory.md, section 6.4). For each
attaching species, per unit area of a face:

- **Deposition.** Cells arrive at the transfer velocity k of
  :mod:`marse.spatial.colloids`, and a fraction alpha of them binds:
  J = alpha k c B(theta). The fraction alpha is the attachment efficiency of
  the species on that face's material, under its conditioning film.
- **Blocking.** Bound cells cover the area they block,
  theta = sum_s n_s a_s, and binding stops as theta approaches the jamming
  limit of random sequential adsorption, theta_J = 0.547 (Feder 1980):
  B = max(0, 1 - theta / theta_J). This is the blocked-area model that
  flow-chamber experiments are fitted with (Busscher and van der Mei 2006).
- **Reversible, then locked.** A bound cell sits in the species' reversible
  component R at first. It detaches back into the liquid at k_off, or locks at
  k_lock into the species' biomass L, as adhesion forces strengthen over the
  first minutes (Mei et al. 2009). Locked cells grow through the network's
  processes; reversible cells last minutes and do not grow.

Deposition and detachment cross the substratum face, so the ledger books them
as imports and exports, as it books the top face. Locking converts R into L,
two components with the same formula, so it conserves everything exactly; the
engine runs it as an extra row of the stoichiometric matrix.

The state is a concentration in the bottom layer of voxels, in mol per m3
(amol per um3). An areal density is that concentration times the height of a
voxel, and a number of cells is an amount divided by the carbon of one cell.
Everything is evaluated at concentrations clamped at zero, as the rates are,
and differentiated consistently, so below zero nothing responds.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

__all__ = ["JAMMING_COVERAGE", "AttachingSpecies", "SurfaceExchange"]

JAMMING_COVERAGE = 0.547
"""The jamming limit of random sequential adsorption of disks (Feder 1980)."""

type Field = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class AttachingSpecies:
    """One species that binds to the substratum, with its material-dependent rates per face."""

    reversible: int  # component index of reversibly bound cells
    attached: int  # component index of locked cells: the species' biomass
    arrival_um_per_h: float  # the transfer velocity k
    cells_per_um3: float  # the suspension in the bulk liquid
    amol_per_cell: float  # carbon per cell, in the state's units
    blocked_area_um2: float
    efficiency: Field  # alpha, per face
    detachment_per_h: Field  # k_off, per face
    locking_per_h: Field  # k_lock, per face


class SurfaceExchange:
    """Binding, detachment and locking of every attaching species at the substratum."""

    def __init__(
        self, species: Sequence[AttachingSpecies], components: int, height_um: float
    ) -> None:
        self.species = tuple(species)
        self.components = components
        self.height_um = height_um
        rows = np.zeros((len(self.species), components))
        for s, sp in enumerate(self.species):
            rows[s, sp.reversible] = -1.0
            rows[s, sp.attached] = 1.0
        self.locking_rows = rows
        # Binding per face on a bare surface, in amol per um2 per h, per species.
        self._binding = [
            sp.efficiency * sp.arrival_um_per_h * sp.cells_per_um3 * sp.amol_per_cell
            for sp in self.species
        ]

    # -- what the surface holds ------------------------------------------------------------

    def cells_per_um2(self, c: Field, s: int) -> Field:
        """Species s bound to each face, reversibly or locked, in cells per um2."""
        sp = self.species[s]
        bottom = np.maximum(c[[sp.reversible, sp.attached], ..., 0], 0.0).sum(axis=0)
        return bottom * self.height_um / sp.amol_per_cell

    def coverage(self, c: Field) -> Field:
        """The area fraction blocked by bound cells, theta, on each face."""
        theta = np.zeros(c.shape[1:-1])
        for s, sp in enumerate(self.species):
            theta = theta + self.cells_per_um2(c, s) * sp.blocked_area_um2
        return theta

    # -- the exchange across the substratum face -------------------------------------------

    def exchange(self, c: Field) -> Field:
        """Into each bottom voxel through its substratum face, per component: mol m-3 um h-1.

        The units are those of a face flux, so that dividing by the voxel
        height gives a rate of change, and summing over faces an import.
        """
        free = np.maximum(1.0 - self.coverage(c) / JAMMING_COVERAGE, 0.0)
        flux = np.zeros((self.components, *c.shape[1:-1]))
        for s, sp in enumerate(self.species):
            reversible = np.maximum(c[sp.reversible, ..., 0], 0.0)
            detached = sp.detachment_per_h * reversible * self.height_um
            flux[sp.reversible] += self._binding[s] * free - detached
        return flux

    def exchange_jacobian(self, c: Field) -> Field:
        """d(exchange / height) / dc in the bottom voxel: shape (J, J, *faces)."""
        lateral = c.shape[1:-1]
        jacobian = np.zeros((self.components, self.components, *lateral))
        blocking = self.coverage(c) < JAMMING_COVERAGE  # past jamming, binding has no slope
        for s, sp in enumerate(self.species):
            for other in self.species:
                # Blocking: binding falls as either state of any species covers the face.
                slope = self._binding[s] / JAMMING_COVERAGE * other.blocked_area_um2
                slope /= other.amol_per_cell
                for k in (other.reversible, other.attached):
                    live = c[k, ..., 0] >= 0
                    jacobian[sp.reversible, k] -= np.where(blocking & live, slope, 0.0)
            live = c[sp.reversible, ..., 0] >= 0
            jacobian[sp.reversible, sp.reversible] -= np.where(live, sp.detachment_per_h, 0.0)
        return jacobian

    # -- locking: an extra process per species, in the bottom layer -----------------------

    def locking(self, c: Field) -> Field:
        """The locking rate of every species in every voxel: shape (S, *voxels)."""
        rates = np.zeros((len(self.species), *c.shape[1:]))
        for s, sp in enumerate(self.species):
            rates[s, ..., 0] = sp.locking_per_h * np.maximum(c[sp.reversible, ..., 0], 0.0)
        return rates

    def locking_jacobian(self, c: Field) -> Field:
        """d(locking) / dc: shape (S, J, *voxels)."""
        jacobian = np.zeros((len(self.species), self.components, *c.shape[1:]))
        for s, sp in enumerate(self.species):
            live = c[sp.reversible, ..., 0] >= 0
            jacobian[s, sp.reversible, ..., 0] = np.where(live, sp.locking_per_h, 0.0)
        return jacobian
