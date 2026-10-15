"""The box's top face open to the air: gases dissolve and escape across it.

Where the top of the box is at the air, as the surface of a salivary film or
of a colony biofilm is, each gas the domain lists crosses the top face and
nothing else does. The air holds the face at the gas's saturation C_s, its
concentration in water in equilibrium with the air (docs/theory.md, sections
4.4 and 4.10). Across the half voxel between the top voxel's centre and the
face, the flux into the box is

    2 D / h (C_s - c_top)

per unit area, for a voxel of height h. That is the transfer a face of fixed
concentration has in the finite-volume scheme (section 4.7), and it runs as
two processes in the top voxel, so that the limiter and the ledger treat it as
they treat any other: the air gives at k C_s and takes back at k c_top, with
k = 2 D / h^2. What the two exchange is an import, booked as exchanged with
the air.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

__all__ = ["AirExchange"]

type Field = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class AirExchange:
    """Gases exchanged with the air through the top face of a box.

    ``gases`` are component indices, each held at ``saturation_mol_per_m3`` at
    the face; ``rate_per_h`` is k = 2 D / h^2 for each. ``shape`` is the box's
    and ``components`` the network's count.
    """

    gases: NDArray[np.intp]
    saturation_mol_per_m3: NDArray[np.float64]
    rate_per_h: NDArray[np.float64]
    shape: tuple[int, ...]
    components: int
    _top: Field = field(init=False, repr=False)

    def __post_init__(self) -> None:
        top = np.zeros(self.shape)
        top[..., -1] = 1.0
        object.__setattr__(self, "_top", top)

    @classmethod
    def at(
        cls,
        gases: NDArray[np.intp],
        saturation_mol_per_m3: NDArray[np.float64],
        diffusivity_um2_per_h: NDArray[np.float64],
        voxel_um: float,
        shape: tuple[int, ...],
    ) -> AirExchange:
        """The exchange for a box of ``shape`` in voxels of ``voxel_um``, from each diffusivity."""
        return cls(
            gases,
            saturation_mol_per_m3,
            2.0 * diffusivity_um2_per_h[gases] / voxel_um**2,
            shape,
            diffusivity_um2_per_h.size,
        )

    @property
    def rows(self) -> NDArray[np.float64]:
        """Two processes per gas, as rows of the stoichiometric matrix: from the air, back to it."""
        count = self.gases.size
        rows = np.zeros((2 * count, self.components))
        rows[np.arange(count), self.gases] = 1.0
        rows[count + np.arange(count), self.gases] = -1.0
        return rows

    def rates(self, c: Field) -> Field:
        """Each process's rate in every voxel: k C_s and k c in the top layer, zero below."""
        dims = len(self.shape)
        k = self.rate_per_h.reshape((-1,) + (1,) * dims)
        given = k * self.saturation_mol_per_m3.reshape((-1,) + (1,) * dims) * self._top
        taken = k * np.maximum(c[self.gases], 0.0) * self._top
        return np.concatenate((given, taken))

    def jacobian(self, c: Field) -> Field:
        """d rate / d c, shape (processes, components, *voxels)."""
        count = self.gases.size
        jacobian = np.zeros((2 * count, self.components, *self.shape))
        live = c[self.gases] >= 0
        dims = len(self.shape)
        k = self.rate_per_h.reshape((-1,) + (1,) * dims)
        jacobian[count + np.arange(count), self.gases] = np.where(live, k * self._top, 0.0)
        return jacobian

    def exchanged(self, extents: Field) -> NDArray[np.float64]:
        """What the air gave less what it took, per component, summed over the box."""
        count = self.gases.size
        net = (extents[:count] - extents[count:]).reshape(count, -1).sum(axis=1)
        exchanged = np.zeros(self.components)
        exchanged[self.gases] = net
        return exchanged
