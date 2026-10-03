"""Spreading: biomass that outgrows its voxel pushes the excess into its neighbours.

Every particulate component that takes up room has a packing density rho_j,
its concentration when it alone fills a voxel. A voxel's biomass fills the
fraction

    phi = sum_j c_j / rho_j,

and no voxel may hold more than it can, phi <= 1, every species counted
together (docs/theory.md, section 6.1). Growth runs in place for a spreading
interval; spreading then moves the excess on, as the continuum model of
Alpkvist and Klapper (2007) does, written on voxels (docs/theory.md, section
9.10):

1. **The full region** Omega holds the voxels with phi >= 1, and e = phi - 1
   is each one's excess.
2. **Pressure.** The discrete Poisson problem sum_b (p_a - p_b) = e_a is
   solved on Omega, with p = 0 in every voxel that has room, and no flux
   through the substratum or the top face. The volume crossing the face from
   a to b is q_ab = p_a - p_b, in voxel volumes.
3. **The sweep.** Voxels are visited in order of decreasing pressure. Each
   has by then received everything flowing into it, and sends q_ab to each
   lower neighbour, carrying every component that takes up room in
   proportion to what it now holds.
4. **Rounds.** A voxel with room may receive more than its room. It then joins
   Omega, and 1 to 3 repeat until every voxel fits.

Every transfer leaves one voxel and enters another, so the sweep conserves
each component exactly, and a voxel never sends more than it holds, so it is
positive at any interval. Components that the spreading carries (such as the
buffer of cell walls) move with the biomass without taking room. Other
components without a density stay where they are, and dissolved components
are left to diffusion.

This increment (Stage 2d, 2d.2) solves the pressure in a column. Boxes in two
and three dimensions come with increment 2d.3.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "CAPACITY_TOLERANCE",
    "ContinuumSpreading",
    "SpreadStats",
    "SpreadingError",
    "volume_fraction",
]

CAPACITY_TOLERANCE = 1e-12
"""How far above 1 a voxel's volume fraction may sit, from rounding alone."""


class SpreadingError(RuntimeError):
    """The biomass cannot be spread: the box is full, or it has reached the top layer."""


def volume_fraction(state: NDArray[np.float64], densities: NDArray[np.float64]) -> NDArray:
    """The fraction of each voxel the biomass fills, sum_j c_j / rho_j.

    ``state`` has shape (components, *voxels); ``densities`` holds rho_j in mol
    per m3, or 0 for a component that takes no room.
    """
    room = densities > 0
    inverse = np.zeros_like(densities, dtype=float)
    inverse[room] = 1.0 / densities[room]
    return np.tensordot(inverse, np.asarray(state, dtype=float), axes=(0, 0))


@dataclass(frozen=True, slots=True)
class SpreadStats:
    """What one spreading step did."""

    rounds: int
    moved: float  # voxel volumes carried across faces, summed
    largest_before: float  # the fullest voxel's phi before spreading


class ContinuumSpreading:
    """Spreading by the pressure of the excess: Darcy flow on voxels, swept downhill."""

    name = "continuum_pressure"
    version = "continuum_pressure_v1"

    def __init__(
        self,
        shape: tuple[int, ...],
        densities: NDArray[np.float64],
        carried: NDArray[np.bool_] | None = None,
    ) -> None:
        if len(shape) != 1:
            raise ValueError("spreading in two and three dimensions arrives in increment 2d.3")
        self.shape = shape
        self.densities = np.asarray(densities, dtype=float)
        # What moves: what takes room, and what it carries along without taking any.
        self.moving = self.densities > 0
        if carried is not None:
            self.moving = self.moving | np.asarray(carried, dtype=bool)

    def spread(self, state: NDArray[np.float64]) -> tuple[NDArray[np.float64], SpreadStats]:
        """Move every voxel's excess on until each fits; returns the new state."""
        state = np.array(state, dtype=float)
        phi = volume_fraction(state, self.densities)
        largest = float(phi.max()) if phi.size else 0.0
        rounds, moved = 0, 0.0
        while phi.max() > 1.0 + CAPACITY_TOLERANCE:
            rounds += 1
            if rounds > phi.size + 1:
                raise SpreadingError("spreading did not settle; this is a bug in MARSE")
            full = phi >= 1.0 - CAPACITY_TOLERANCE
            if full.all():
                raise SpreadingError(
                    "the box is full: every voxel holds all the biomass it can, and there is "
                    "nowhere for more to go"
                )
            pressure = _column_pressure(full, np.maximum(phi - 1.0, 0.0))
            moved += self._sweep(state, phi, pressure)
            phi = volume_fraction(state, self.densities)
        return state, SpreadStats(rounds, moved, largest)

    def _sweep(
        self, state: NDArray[np.float64], phi: NDArray[np.float64], pressure: NDArray[np.float64]
    ) -> float:
        """Send each face's flow downhill, voxel by voxel in order of falling pressure.

        Material keeps its order. A voxel holds a stack of parcels, each a
        volume and the amounts in it: its own content, with what flows in from
        below put beneath it and what flows in from above put on top. What
        leaves through the top face is taken from the top of the stack, and
        what leaves through the bottom face from its bottom. A layer of cells
        therefore stays a layer as it moves, mixed only within the voxels it
        straddles, as Wanner and Gujer's displacement carries it. Mixing every
        voxel before it sends instead (donor cell) smeared a labelled band over
        half the film in three doublings (docs/validation.md).
        """
        n = pressure.size
        flow = pressure[:-1] - pressure[1:]  # across the face above voxel i, upwards positive
        moving = state[self.moving]
        stacks: list[list[_Parcel]] = [
            [_Parcel(float(phi[v]), moving[:, v].copy())] for v in range(n)
        ]
        sent_total = 0.0
        # A stable sort, so voxels of equal pressure, which exchange nothing, keep their order.
        for v in map(int, np.argsort(-pressure, kind="stable")):
            if pressure[v] <= 0.0:
                break  # the rest have room and only receive
            if v + 1 < n and flow[v] > 0.0:
                parcels = _take(stacks[v], float(flow[v]), from_top=True)
                stacks[v + 1][:0] = parcels  # enters the voxel above through its floor
                sent_total += float(flow[v])
            if v > 0 and flow[v - 1] < 0.0:
                parcels = _take(stacks[v], float(-flow[v - 1]), from_top=False)
                stacks[v - 1].extend(parcels)  # enters the voxel below through its ceiling
                sent_total += float(-flow[v - 1])
        for v, stack in enumerate(stacks):
            total = np.zeros(moving.shape[0])
            for parcel in stack:
                total = total + parcel.amounts
            moving[:, v] = total
        state[self.moving] = moving
        return sent_total


@dataclass(slots=True)
class _Parcel:
    """A slice of material in a voxel: its volume, in voxel volumes, and the amounts in it."""

    volume: float
    amounts: NDArray[np.float64]


def _take(stack: list[_Parcel], volume: float, *, from_top: bool) -> list[_Parcel]:
    """Remove ``volume`` from one end of the stack; returns it, bottom first.

    A parcel cut in two keeps the amounts in proportion to volume; the part
    left behind is what remains after subtraction, so nothing is lost, and
    since the part taken is at most the whole, nothing goes negative.
    """
    taken: list[_Parcel] = []
    wanted = volume
    while wanted > 0.0 and stack:
        parcel = stack[-1] if from_top else stack[0]
        if parcel.volume <= wanted:
            stack.pop(-1 if from_top else 0)
            taken.append(parcel)
            wanted -= parcel.volume
            continue
        fraction = wanted / parcel.volume
        cut = _Parcel(wanted, parcel.amounts * fraction)
        parcel.amounts = parcel.amounts - cut.amounts
        parcel.volume -= wanted
        taken.append(cut)
        wanted = 0.0
    if from_top:
        taken.reverse()  # taken top first; returned bottom first
    return taken


def _column_pressure(full: NDArray[np.bool_], excess: NDArray[np.float64]) -> NDArray:
    """The pressure in a column: sum_b (p_a - p_b) = e_a on the full voxels, 0 elsewhere.

    Neighbours are the voxels above and below; the substratum and the top face
    pass nothing, so the end voxels have one neighbour. The tridiagonal system
    is solved by the Thomas algorithm, with an identity row for every voxel
    that has room.
    """
    n = full.size
    lower = np.zeros(n)
    diagonal = np.ones(n)
    upper = np.zeros(n)
    rhs = np.zeros(n)
    for i in map(int, np.flatnonzero(full)):
        neighbours = int(i > 0) + int(i < n - 1)
        diagonal[i] = neighbours
        if i > 0 and full[i - 1]:
            lower[i] = -1.0
        if i < n - 1 and full[i + 1]:
            upper[i] = -1.0
        rhs[i] = excess[i]
    c, d = np.zeros(n), np.zeros(n)
    c[0], d[0] = upper[0] / diagonal[0], rhs[0] / diagonal[0]
    for i in range(1, n):
        beta = diagonal[i] - lower[i] * c[i - 1]
        c[i] = upper[i] / beta
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / beta
    pressure = np.empty(n)
    pressure[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        pressure[i] = d[i] - c[i] * pressure[i + 1]
    return pressure
