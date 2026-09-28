"""A well-mixed reservoir bordering the box: the mouth's saliva over a site of plaque.

A salivary film, the top voxels of the box, is renewed from a pool of liquid
at a rate k in each of its voxels, and the pool's own composition changes
with what the film returns, what enters the pool, and what leaves it
(docs/theory.md, sections 4.8 and 9.9). The box's top face is closed: the
film's surface is open to the air.

The S1 prototype first solved the pool apart from the box, and corrected it
after each span by exactly what had crossed. That conserved to rounding, but
while an intake mixed the film fast, the answer changed by 0.03 pH between
spans of 60 s and 5 s. So the pool is solved with the box, implicitly:

- **Its unknowns** are u = n / H_ref: the pool's amount of each component per
  unit area of substratum, over a fixed reference thickness. They look like
  concentrations to the error control and to Newton's method, and they are
  amounts to the update, which books exactly what crosses.
- **Its thickness** H(t), the pool's volume per unit area of substratum, is
  given over each span (:class:`ReservoirPath`), with what enters it. The
  pool's concentration is m = u H_ref / H(t).
- **The exchange** into voxel v is k_v (m - c_v) for each exchanged component.
  It runs as two processes per component: one brings the pool's liquid in,
  one takes the voxel's back. The limiter scales what leaves a voxel like any
  consumption. What the pool gives is checked after, and scaled down if the
  pool would go below zero.
- **The linear systems** are the box's, bordered by the pool:
  [S, -aB; -aC, 1 - aD]. The Schur complement on the pool needs the box's
  response to each exchanged component of the pool, one solve each per
  matrix. Every solve after that costs one box solve and one J x J product.

The state holds both, shape ``(J, voxels + 1)``: the box flattened, then the
pool. :meth:`ReservoirTransport.pack` and :meth:`~ReservoirTransport.unpack`
convert.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse.core.implicit import (  # the engine's own step, reused
    GAMMA,
    ReactionTransport,
    StepStats,
    Tolerance,
    _magnitude,
    _Stage,
)
from marse.microbes.adhesion import SurfaceExchange
from marse.spatial.column import ColumnSystem
from marse.spatial.multigrid import ImplicitSystem
from marse.spatial.transport import Diffusion, divergence

__all__ = ["BorderedSystem", "ReservoirPath", "ReservoirTransport"]

type Field = NDArray[np.float64]

_MARGIN = 1.0 - 1e-12
_RESPONSE_TOLERANCE = 1e-10  # the box's response to the pool, where it is solved iteratively


@dataclass(frozen=True, slots=True)
class ReservoirPath:
    """The pool's thickness over one span, and what enters it.

    ``thickness_um`` and ``growth_um_per_h`` sample H and dH/dt at
    ``times_h``, counted from ``start_h``; between samples H is the cubic
    Hermite polynomial through them, so it is smooth inside the span. A span
    ends at every event that changes H abruptly, such as a swallow.

    Liquid enters with the growth of H: ``drink_um_per_h`` of it is a drink of
    composition ``drink_mol_per_m3``, the rest is secretion of composition
    ``secreted_mol_per_m3``. ``supply_per_h`` adds amounts without volume, such
    as sugar dissolving from a sweet, per unit area per hour, in mol/m3 x um.
    """

    start_h: float
    times_h: NDArray[np.float64]
    thickness_um: NDArray[np.float64]
    growth_um_per_h: NDArray[np.float64]
    secreted_mol_per_m3: NDArray[np.float64]
    drink_um_per_h: float = 0.0
    drink_mol_per_m3: NDArray[np.float64] | None = None
    supply_per_h: NDArray[np.float64] | None = None

    def at(self, t_h: float) -> tuple[float, float]:
        """H and dH/dt at time ``t_h``."""
        times = self.times_h
        if times.size == 1:
            return float(self.thickness_um[0]), float(self.growth_um_per_h[0])
        x = t_h - self.start_h
        i = int(np.clip(np.searchsorted(times, x) - 1, 0, times.size - 2))
        width = times[i + 1] - times[i]
        u = min(max((x - times[i]) / width, 0.0), 1.0)
        h0, h1 = self.thickness_um[i], self.thickness_um[i + 1]
        d0, d1 = self.growth_um_per_h[i] * width, self.growth_um_per_h[i + 1] * width
        value = (
            (2 * u**3 - 3 * u**2 + 1) * h0
            + (u**3 - 2 * u**2 + u) * d0
            + (-2 * u**3 + 3 * u**2) * h1
            + (u**3 - u**2) * d1
        )
        slope = (
            (6 * u**2 - 6 * u) * h0
            + (3 * u**2 - 4 * u + 1) * d0
            + (-6 * u**2 + 6 * u) * h1
            + (3 * u**2 - 2 * u) * d1
        ) / width
        return float(value), float(slope)

    def _sources(self, grown: float, span: float) -> NDArray[np.float64]:
        drunk = self.drink_um_per_h * span
        added = self.secreted_mol_per_m3 * (grown - drunk)
        if self.drink_mol_per_m3 is not None:
            added = added + self.drink_mol_per_m3 * drunk
        if self.supply_per_h is not None:
            added = added + self.supply_per_h * span
        return added

    def rate(self, t_h: float) -> NDArray[np.float64]:
        """What enters per unit area per hour at ``t_h``, per component, in mol/m3 x um/h."""
        return self._sources(self.at(t_h)[1], 1.0)

    def added(self, t0_h: float, t1_h: float) -> NDArray[np.float64]:
        """What enters per unit area between two times, exactly as the thickness grows."""
        return self._sources(self.at(t1_h)[0] - self.at(t0_h)[0], t1_h - t0_h)


class BorderedSystem:
    """The box's implicit system bordered by the pool, solved by the Schur complement.

    With S the box's matrix, B how the box's rates respond to the pool, C how
    the pool's rates respond to the box and D to itself, all per component,

        S x_b - a B x_u = r_b,     -a C x_b + (1 - a D) x_u = r_u,

    gives x_b = S^-1 r_b + a W x_u with W = S^-1 B, and a J x J system for x_u.
    """

    def __init__(
        self,
        box: ColumnSystem | ImplicitSystem,
        a: float,
        exchanged: NDArray[np.intp],
        to_box: Field,
        to_pool: Field,
        own: NDArray[np.float64],
        shape: tuple[int, ...],
    ) -> None:
        self.box, self.a, self.exchanged, self.shape = box, a, exchanged, shape
        self._to_box, self._own = to_box, own
        count = exchanged.size
        j = int(box.diffusivity.size)
        self.components = j
        self.voxels = math.prod(shape)
        self.to_pool = to_pool.reshape(count, -1)
        responses = np.zeros((count, j, self.voxels))
        for i, component in enumerate(exchanged):
            unit = np.zeros((j, *shape))
            unit[component] = to_box[i]
            responses[i] = self._box_solve(unit).reshape(j, -1)
        self.responses = responses
        schur = np.eye(j)
        schur[exchanged, exchanged] = 1.0 - a * own
        coupling = np.einsum("iv,kiv->ik", self.to_pool, responses[:, exchanged, :])
        schur[np.ix_(exchanged, exchanged)] -= a * a * coupling
        self.schur_inverse = np.linalg.inv(schur)

    def apply(self, x: Field) -> Field:
        """The bordered matrix applied to a packed state: what :meth:`solve` inverts."""
        j = self.components
        box = x[:, : self.voxels].reshape(j, *self.shape)
        pool = x[:, self.voxels]
        applied = self.box.apply(box).reshape(j, -1)
        count = self.exchanged.size
        applied[self.exchanged] -= (
            self.a * self._to_box.reshape(count, -1) * pool[self.exchanged, None]
        )
        pool_out = pool.copy()
        pool_out[self.exchanged] = (1.0 - self.a * self._own) * pool[self.exchanged]
        pool_out[self.exchanged] -= self.a * np.einsum(
            "iv,iv->i", self.to_pool, box[self.exchanged].reshape(count, -1)
        )
        return np.concatenate((applied, pool_out[:, None]), axis=1)

    def _box_solve(self, b: Field) -> Field:
        if isinstance(self.box, ColumnSystem):
            return self.box.precondition(b)
        scale = np.full(b.shape, max(float(np.abs(b).max()), 1e-300))
        solved, _ = self.box.solve(b, scale=scale, tolerance=_RESPONSE_TOLERANCE)
        return solved

    def solve(
        self,
        b: Field,
        *,
        scale: Field,
        tolerance: float = 1e-6,
        max_iterations: int = 200,
    ) -> tuple[Field, int]:
        """Solve for the box and the pool together; ``b`` and the result are packed states."""
        j, count = self.components, self.exchanged.size
        r_box = b[:, : self.voxels].reshape(j, *self.shape)
        box_scale = np.broadcast_to(scale, b.shape)[:, : self.voxels].reshape(j, *self.shape)
        z, iterations = self.box.solve(
            r_box, scale=box_scale, tolerance=tolerance, max_iterations=max_iterations
        )
        rhs = b[:, self.voxels].copy()
        rhs[self.exchanged] += self.a * np.einsum(
            "iv,iv->i", self.to_pool, z[self.exchanged].reshape(count, -1)
        )
        pool = self.schur_inverse @ rhs
        box = z.reshape(j, -1) + self.a * np.einsum(
            "i,ijv->jv", pool[self.exchanged], self.responses
        )
        return np.concatenate((box, pool[:, None]), axis=1), iterations


class ReservoirTransport(ReactionTransport):
    """A reaction network in a box closed at the top, bordered by a well-mixed pool.

    ``exchanged`` lists the components the pool and the box exchange: the
    dissolved ones. ``reference_um`` is H_ref, the thickness that turns the
    pool's amounts into its unknowns. Before each span, the caller sets
    :attr:`path` and :attr:`exchange_per_h`, the rate k in every voxel.
    """

    def __init__(
        self,
        diffusion: Diffusion,
        stoichiometry: NDArray[np.float64],
        rates: Callable[[Field], Field],
        jacobian: Callable[[Field], Field],
        *,
        exchanged: Sequence[int],
        reference_um: float,
        surface: SurfaceExchange | None = None,
    ) -> None:
        if not diffusion.closed_top:
            raise ValueError("a box under a film exchanges through the film: close its top face")
        super().__init__(diffusion, stoichiometry, rates, jacobian, surface)
        self.exchanged = np.asarray(exchanged, dtype=np.intp)
        count = self.exchanged.size
        components = self.stoichiometry.shape[1]
        self._reacting = self.stoichiometry.shape[0]
        rows = np.zeros((2 * count, components))
        rows[np.arange(count), self.exchanged] = 1.0  # from the pool into a voxel
        rows[count + np.arange(count), self.exchanged] = -1.0  # from a voxel back to the pool
        self.stoichiometry = np.vstack((self.stoichiometry, rows))
        self.consumed = np.maximum(-self.stoichiometry, 0.0)
        self.produced = np.maximum(self.stoichiometry, 0.0)
        self.reference_um = float(reference_um)
        grid = diffusion.grid
        self.voxels = grid.voxels
        self.areal_um = grid.voxel_volume_um3 / grid.footprint_um2
        self.exchange_per_h = np.zeros(self.shape)
        self.path: ReservoirPath | None = None
        self._time = 0.0

    # -- the packed state ----------------------------------------------------------------

    def pack(self, box: Field, pool: NDArray[np.float64]) -> Field:
        return np.concatenate((box.reshape(box.shape[0], -1), pool[:, None]), axis=1)

    def unpack(self, y: Field) -> tuple[Field, NDArray[np.float64]]:
        return y[:, : self.voxels].reshape(y.shape[0], *self.shape), y[:, self.voxels]

    def _per_area(self, field: Field) -> NDArray[np.float64]:
        """Summed over the box, per unit area of substratum: mol/m3 x um."""
        return field.reshape(field.shape[0], -1).sum(axis=1) * self.areal_um

    def _pool_concentration(self, pool: NDArray[np.float64]) -> tuple[NDArray, float]:
        assert self.path is not None, "set the path before each span"
        thickness, _ = self.path.at(self._time)
        return pool * self.reference_um / thickness, thickness

    # -- the right-hand side and its Jacobian --------------------------------------------

    def evaluate(self, y: Field) -> _Stage:
        box, pool = self.unpack(y)
        concentration, _ = self._pool_concentration(pool)
        fluxes = self.diffusion.fluxes(box)
        k = self.exchange_per_h
        dims = len(self.shape)
        given = np.maximum(concentration[self.exchanged], 0.0).reshape((-1,) + (1,) * dims)
        inflow = k * given
        outflow = k * np.maximum(box[self.exchanged], 0.0)
        reactions = np.concatenate((self.rates(box), inflow, outflow))
        rate = divergence(fluxes, self.spacing) + self._react(reactions)
        substratum = None
        if self.surface is not None:
            substratum = self.surface.exchange(box)
            rate[..., 0] += substratum / self.spacing[-1]
        assert self.path is not None
        pool_rate = self.path.rate(self._time) / self.reference_um
        pool_rate[self.exchanged] -= self._per_area(inflow - outflow) / self.reference_um
        return _Stage(self.pack(rate, pool_rate), fluxes, reactions, substratum)

    def _system(self, y: Field, a: float) -> BorderedSystem:  # type: ignore[override]
        box, pool = self.unpack(y)
        concentration, thickness = self._pool_concentration(pool)
        count = self.exchanged.size
        jacobian = self.jacobian(box)
        extra = np.zeros((2 * count, *jacobian.shape[1:]))
        live = box[self.exchanged] >= 0
        extra[count + np.arange(count), self.exchanged] = np.where(live, self.exchange_per_h, 0.0)
        blocks = np.einsum("pj,pk...->jk...", self.stoichiometry, np.concatenate((jacobian, extra)))
        if self.surface is not None:
            blocks[..., 0] += self.surface.exchange_jacobian(box)
        solver = ColumnSystem if len(self.shape) == 1 else ImplicitSystem
        system = solver(
            self.shape,
            self.spacing,
            self.diffusion.diffusivity_um2_per_h,
            a,
            blocks,
            closed_top=True,
        )
        k = self.exchange_per_h
        alive = concentration[self.exchanged] >= 0
        dims = len(self.shape)
        to_box = np.where(
            alive.reshape((-1,) + (1,) * dims), k * self.reference_um / thickness, 0.0
        )
        to_pool = np.where(live, k * self.areal_um / self.reference_um, 0.0)
        own = np.where(alive, -float(k.sum()) * self.areal_um / thickness, 0.0)
        return BorderedSystem(system, a, self.exchanged, to_box, to_pool, own, self.shape)

    # -- one step --------------------------------------------------------------------------

    def _conserve(
        self,
        box: Field,
        pool: NDArray[np.float64],
        added: NDArray[np.float64],
        transfers: tuple[Field, ...],
        extents: Field,
        substratum: Field | None,
    ) -> tuple[Field, NDArray[np.float64], Field, Field | None, int]:
        """The flux-form update of the box and the pool, both kept non-negative.

        The box is limited as the engine limits it. The pool then gains what
        entered it and what the film returned, and loses what it gave; if it
        would fall below zero in a component, what it gives of that component
        is scaled down, and the box is updated again.
        """
        count, first = self.exchanged.size, self._reacting
        dims = len(self.shape)
        for _ in range(51):
            new_box, transfers, extents, substratum, rounds = self._limited_update(
                box, transfers, extents, substratum
            )
            given = self._per_area(extents[first : first + count])
            returned = self._per_area(extents[first + count :])
            new_pool = pool + added / self.reference_um
            new_pool[self.exchanged] += (returned - given) / self.reference_um
            short = new_pool[self.exchanged] < 0
            if not short.any():
                return new_box, new_pool, extents, substratum, rounds
            available = pool[self.exchanged] * self.reference_um + added[self.exchanged] + returned
            share = np.where(
                short, np.maximum(available, 0.0) * _MARGIN / np.where(short, given, 1.0), 1.0
            )
            extents = extents.copy()
            extents[first : first + count] *= share.reshape((-1,) + (1,) * dims)
        raise ArithmeticError("the pool could not be kept from going negative")  # pragma: no cover

    def _step(
        self, y: Field, h: float, reference: Field, atol: Tolerance, rtol: float
    ) -> tuple[Field, NDArray[np.float64], Field, int, int]:
        g = GAMMA
        start = self._now
        reference = _magnitude(reference)
        self._time = start
        system = self._system(y, g * h)
        self._time = start + g * h
        stage1, system, n1 = self._stage(y, g * h, y, system, reference, atol, rtol)
        s1 = self.evaluate(stage1)
        base = y + (1 - g) * h * s1.rate
        self._time = start + h
        stage2, system, n2 = self._stage(base, g * h, stage1, system, reference, atol, rtol)
        s2 = self.evaluate(stage2)
        transfers = tuple(
            h * ((1 - g) * f1 + g * f2) for f1, f2 in zip(s1.fluxes, s2.fluxes, strict=True)
        )
        extents = h * ((1 - g) * s1.reactions + g * s2.reactions)
        substratum = None
        if s1.substratum is not None and s2.substratum is not None:
            substratum = h * ((1 - g) * s1.substratum + g * s2.substratum)
        box, pool = self.unpack(y)
        assert self.path is not None
        added = self.path.added(start, start + h)
        new_box, new_pool, extents, substratum, rounds = self._conserve(
            box, pool, added, transfers, extents, substratum
        )
        new = self.pack(new_box, new_pool)
        scale = atol + rtol * np.maximum(reference, _magnitude(y, new))
        estimate, _ = system.solve(g * h * (s2.rate - s1.rate), scale=scale, tolerance=1e-2)
        count, first = self.exchanged.size, self._reacting
        exchanged = extents[first : first + count] - extents[first + count :]
        imports = np.zeros(y.shape[0])
        imports[self.exchanged] = exchanged.reshape(count, -1).sum(axis=1)
        if substratum is not None:  # bound less detached, as in the engine
            imports += substratum.reshape(y.shape[0], -1).sum(axis=1) / self.spacing[-1]
        return new, imports, estimate, rounds, n1 + n2

    def integrate(  # type: ignore[override]
        self,
        y: Field,
        span: float,
        *,
        relative_tolerance: float,
        absolute_tolerance: Tolerance,
        first_step: float | None = None,
        peak: Field | None = None,
        start_h: float = 0.0,
    ) -> tuple[Field, NDArray[np.float64], StepStats]:
        """Advance the packed state by ``span`` from ``start_h``; imports are the box's."""
        self._time = start_h
        return super().integrate(
            y,
            span,
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
            first_step=first_step,
            peak=peak,
            start_h=start_h,
        )
