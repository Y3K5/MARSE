"""Implicit, conservative, positive integration of reactions and transport together.

On a grid, concentrations change as

    dc/dt = divergence of the face fluxes + N^T r(c),

diffusion (:mod:`marse.spatial.transport`) and reactions
(:mod:`marse.microbes.kinetics`) at once. Diffusion at physical diffusivities is
stiff: an explicit step on 2 um voxels would have to be a quarter of a
millisecond (docs/theory.md, section 9.1). Splitting reactions from transport
fails differently: at practical steps it starves a biofilm, because each step
supplies one refill of the box instead of a continuous flux. A quasi-steady
treatment fails in exactly the transients MARSE is built to study. So
everything is integrated together, implicitly (docs/theory.md, section 9.8):

- **The scheme** is the two-stage SDIRK method of Alexander (1977), with
  gamma = 1 - 1/sqrt(2):

      Y1 = y + gamma h f(Y1),
      Y2 = y + h ((1 - gamma) f(Y1) + gamma f(Y2)).

  It is L-stable and stiffly accurate, and second order: long steps land on
  the quasi-steady state, and short steps follow a transient.
- **Conservation for any solver tolerance.** The new state is not taken from
  the stage solution. It is rebuilt in flux form, from the face transfers and
  process extents the two stages imply,

      y+ = y + transfers + N^T extents,

  and each transfer leaves one voxel and enters another, while each process
  row conserves carbon, nitrogen and electrons. So the balance holds to
  rounding whatever the Newton or multigrid tolerance.
- **Positivity by limiting.** No scheme above first order is positive by
  itself at every step size. Where the update would leave a value negative,
  the transfers out of that voxel and the processes consuming there are
  scaled down, each as a whole, net of what arrives. A single-pass limiter that
  counts only the stock is the guaranteed fallback. Nothing is clipped. The
  update sums what arrives in each voxel and what leaves it separately, so a
  voxel that loses nothing can only gain, even in rounding. Traces below the
  smallest normal number, where no relative margin survives rounding, stop
  giving instead of being scaled.
- **Newton's method, projected.** A value well above the absolute tolerance
  may fall at most to a tenth of itself per iteration, which keeps Newton out
  of the kink of the rate laws at zero while a large transient is resolved. A
  value already below the tolerance may go negative, because the stage
  solution itself can undershoot there. A Newton iteration that does not
  converge rejects the step, which is retried shorter.
- **Error control.** The difference between the second-order result and a
  first-order one, h gamma (f(Y2) - f(Y1)), filtered through the stage matrix
  so that stiff components do not dominate it, sets the substep. It is
  deterministic, so a run replays bit for bit. A run's first step comes from
  the rates and their change over a trial Euler step (Hairer, Norsett and
  Wanner 1993, section II.4): colonies placed in fresh liquid start far from
  their quasi-steady state, and a step that tried to cross that in one go
  would only fail.

Where cells in the liquid bind to the substratum (:mod:`marse.microbes.adhesion`,
docs/theory.md, section 6.4), the bottom face carries a second exchange beside
the top one. Deposition arrives through it and detachment leaves through it, as
face transfers. The limiter scales them like any other, and the ledger counts
them as imports. Locking turns reversibly bound cells into biomass, as an extra
row of the stoichiometric matrix acting in the bottom layer.

Where the top face is at the air (:mod:`marse.spatial.air`, docs/theory.md,
section 4.10), the gases it holds cross it as two more processes per gas in
the top layer, and the ledger counts what they exchange as an import.

Where a column's plaque spreads (:mod:`marse.biofilm.spreading`, docs/theory.md,
section 6.1), every step is followed by packing the solid from the substratum
up, which detaches whatever passes the plaque's maximum height. Wear is taken
off in two halves, before the step and after it (Strang splitting), so the
step stays second order. What detaches leaves the box, and the ledger counts
it as an export.

The linear systems are solved by :mod:`marse.spatial.multigrid`, or, in a
column, directly by :mod:`marse.spatial.column`. One matrix, with the
Jacobian at the start of the step, serves both stages, and it is rebuilt only
when Newton converges slowly.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse.biofilm.spreading import Spreading
from marse.microbes.adhesion import SurfaceExchange
from marse.spatial.air import AirExchange
from marse.spatial.column import ColumnSystem
from marse.spatial.multigrid import ImplicitSystem
from marse.spatial.transport import Diffusion, divergence

__all__ = ["GAMMA", "NewtonFailure", "ReactionTransport", "StepStats"]

GAMMA = 1.0 - 1.0 / math.sqrt(2.0)
_MARGIN = 1.0 - 1e-12
_TINY = float(np.finfo(float).tiny)  # below it, rounding swallows any relative margin
_NEWTON_TOLERANCE = 1e-3  # corrections, in units of the error tolerance
_NEWTON_ITERATIONS = 20
_LINEAR_TOLERANCE = 1e-4

type Field = NDArray[np.float64]
type Transfers = tuple[Field, ...]
type Tolerance = float | NDArray[np.float64]


def _tolerance(absolute: Tolerance, dims: int) -> Tolerance:
    """An absolute tolerance per component shaped to broadcast against a field.

    One number stays that number, so a run with a single tolerance computes
    exactly what it computed before tolerances could be given per component.
    """
    if np.ndim(absolute) == 0:
        return absolute
    return np.asarray(absolute, dtype=float).reshape((-1,) + (1,) * dims)


def _magnitude(*fields: Field) -> Field:
    """Each component's largest absolute value anywhere in the box, shaped to broadcast.

    Errors are measured against a component's own scale across the box, not
    against the value in each voxel: a trace of lactate seeping into the liquid
    far from the colonies is held to the accuracy lactate needs where it
    matters, not resolved ever more finely where it does not.
    """
    largest = np.max([np.abs(f).reshape(f.shape[0], -1).max(axis=1) for f in fields], axis=0)
    return largest.reshape((-1,) + (1,) * (fields[0].ndim - 1))


class NewtonFailure(ArithmeticError):
    """Newton's method did not converge within its iterations; the step is retried shorter."""


@dataclass(frozen=True, slots=True)
class StepStats:
    """How :meth:`ReactionTransport.integrate` got there."""

    accepted: int
    rejected: int
    limited: int
    newton_failures: int
    newton_iterations: int
    next_step: float


@dataclass(frozen=True, slots=True)
class _Stage:
    rate: Field
    fluxes: Transfers
    reactions: Field
    substratum: Field | None = None  # into the bottom voxels through the substratum


class ReactionTransport:
    """A reaction network in a box of voxels, with the bulk liquid held above it."""

    def __init__(
        self,
        diffusion: Diffusion,
        stoichiometry: NDArray[np.float64],
        rates: Callable[[Field], Field],
        jacobian: Callable[[Field], Field],
        surface: SurfaceExchange | None = None,
        spreading: Spreading | None = None,
        air: AirExchange | None = None,
    ) -> None:
        self.diffusion = diffusion
        self.surface = surface
        self.spreading = spreading
        self.air = air
        stoichiometry = np.asarray(stoichiometry, dtype=float)  # processes x components
        if surface is not None:
            # Locking runs as one more process per species, in the bottom layer.
            stoichiometry = np.vstack((stoichiometry, surface.locking_rows))
            locked_rates, locked_jacobian = rates, jacobian

            def rates(c: Field) -> Field:
                return np.concatenate((locked_rates(c), surface.locking(c)))

            def jacobian(c: Field) -> Field:
                return np.concatenate((locked_jacobian(c), surface.locking_jacobian(c)))

        first = stoichiometry.shape[0]
        if air is not None:
            # The air runs as two more processes per gas, in the top layer.
            stoichiometry = np.vstack((stoichiometry, air.rows))
            aired_rates, aired_jacobian = rates, jacobian

            def rates(c: Field) -> Field:
                return np.concatenate((aired_rates(c), air.rates(c)))

            def jacobian(c: Field) -> Field:
                return np.concatenate((aired_jacobian(c), air.jacobian(c)))

        self._air_rows = slice(first, stoichiometry.shape[0])
        self.stoichiometry = stoichiometry
        self.consumed = np.maximum(-self.stoichiometry, 0.0)
        self.produced = np.maximum(self.stoichiometry, 0.0)
        self.rates = rates
        self.jacobian = jacobian
        self.shape = diffusion.grid.shape
        self.spacing = diffusion.spacing
        self._now = 0.0  # the start of the current step, h
        components = stoichiometry.shape[1]
        self.detached = np.zeros(components)  # solid that left the plaque, summed over steps
        self._detaching = np.zeros(components)  # what the step being tried detached
        self.aired = np.zeros(components)  # what the air gave less what it took, summed
        self._airing = np.zeros(components)  # the same, over the step being tried

    # -- the right-hand side --------------------------------------------------------------

    def _react(self, extents: Field) -> Field:
        return np.einsum("pj,p...->j...", self.stoichiometry, extents)

    def evaluate(self, c: Field) -> _Stage:
        fluxes = self.diffusion.fluxes(c)
        reactions = self.rates(c)
        rate = divergence(fluxes, self.spacing) + self._react(reactions)
        if self.surface is None:
            return _Stage(rate, fluxes, reactions)
        substratum = self.surface.exchange(c)
        rate[..., 0] += substratum / self.spacing[-1]
        return _Stage(rate, fluxes, reactions, substratum)

    def _system(self, c: Field, a: float) -> ImplicitSystem | ColumnSystem:
        blocks = np.einsum("pj,pk...->jk...", self.stoichiometry, self.jacobian(c))
        if self.surface is not None:
            blocks[..., 0] += self.surface.exchange_jacobian(c)
        # A column is solved directly; a box by multigrid (docs/theory.md, section 9.7).
        solver = ColumnSystem if len(self.shape) == 1 else ImplicitSystem
        return solver(
            self.shape,
            self.spacing,
            self.diffusion.diffusivity_um2_per_h,
            a,
            blocks,
            closed_top=self.diffusion.closed_top,
        )

    # -- one stage: Y = base + a f(Y) ----------------------------------------------------

    def _stage(
        self,
        base: Field,
        a: float,
        guess: Field,
        system: ImplicitSystem | ColumnSystem,
        reference: Field,
        atol: Tolerance,
        rtol: float,
    ) -> tuple[Field, ImplicitSystem | ColumnSystem, int]:
        y = guess.copy()
        previous = math.inf
        for iteration in range(1, _NEWTON_ITERATIONS + 1):
            residual = y - base - a * self.evaluate(y).rate
            scale = atol + rtol * np.maximum(reference, _magnitude(y))
            correction, _ = system.solve(-residual, scale=scale, tolerance=_LINEAR_TOLERANCE)
            size = float(np.max(np.abs(correction) / scale))
            if size < _NEWTON_TOLERANCE:
                return y + correction, system, iteration
            if size > 0.5 * previous:  # slow: the Jacobian has gone stale
                system = self._system(y, a)
            previous = size
            trial = y + correction
            y = np.where((y > atol) & (trial < 0.1 * y), 0.1 * y, trial)
        raise NewtonFailure(
            f"Newton's method did not converge in {_NEWTON_ITERATIONS} iterations "
            f"(last correction {size:.3g} of the tolerance)"
        )

    # -- the conservative, positive update ------------------------------------------------

    def _flows(
        self, transfers: Transfers, extents: Field, substratum: Field | None = None
    ) -> tuple[Field, Field]:
        """What arrives in each voxel and what leaves it, per component: both non-negative.

        Arrivals are inflow through the faces and production; departures are
        outflow and consumption. Each transfer departs from one voxel and
        arrives in its neighbour, so the two conserve exactly what the
        divergence does. Kept apart, they make the update positive by
        construction wherever nothing leaves: a voxel that loses nothing can
        only gain, even in rounding.
        """
        into = np.einsum("pj,p...->j...", self.produced, extents)
        out = np.einsum("pj,p...->j...", self.consumed, extents)
        dims = len(self.shape)
        for axis, upper in enumerate(transfers):
            h = self.spacing[axis]
            forward = np.maximum(upper, 0.0) / h  # leaves through the upper face
            backward = np.maximum(-upper, 0.0) / h  # enters through the upper face
            out += forward
            into += backward
            if axis < dims - 1:  # periodic: the upper face of i is the lower face of i + 1
                into += np.roll(forward, 1, axis=axis + 1)
                out += np.roll(backward, 1, axis=axis + 1)
            else:  # the substratum face below the first voxel carries only adhesion
                into[..., 1:] += forward[..., :-1]
                out[..., 1:] += backward[..., :-1]
        if substratum is not None:
            h = self.spacing[-1]
            into[..., 0] += np.maximum(substratum, 0.0) / h  # bound from the suspension
            out[..., 0] += np.maximum(-substratum, 0.0) / h  # detached into the liquid
        return into, out

    def _apply(
        self, y: Field, transfers: Transfers, extents: Field, substratum: Field | None = None
    ) -> tuple[Field, Field, Field]:
        """The new state, with what arrived and what left."""
        into, out = self._flows(transfers, extents, substratum)
        return (y - out) + into, into, out

    def _scale(
        self, transfers: Transfers, extents: Field, share: Field, substratum: Field | None = None
    ) -> tuple[Transfers, Field, Field | None]:
        """Scale what leaves voxel i by share[i]: transfers by their source, processes whole.

        What binds from the suspension comes from the bulk, which is never
        short; what detaches leaves the bottom voxel and is scaled with it.
        """
        dims = len(self.shape)
        scaled = []
        for axis, upper in enumerate(transfers):
            if axis < dims - 1:
                source_below = np.roll(share, -1, axis=axis + 1)
            else:  # the bulk above the top face is never short
                source_below = np.concatenate(
                    (share[..., 1:], np.ones_like(share[..., :1])), axis=-1
                )
            scaled.append(np.where(upper > 0, upper * share, upper * source_below))
        limit = np.where(
            (self.consumed > 0).reshape(self.consumed.shape + (1,) * dims),
            share[None],
            1.0,
        ).min(axis=1)
        if substratum is not None:
            substratum = np.where(substratum < 0, substratum * share[..., 0], substratum)
        return tuple(scaled), extents * limit, substratum

    def _limited_update(
        self, y: Field, transfers: Transfers, extents: Field, substratum: Field | None = None
    ) -> tuple[Field, Transfers, Field, Field | None, int]:
        new, into, out = self._apply(y, transfers, extents, substratum)
        if np.all(new >= 0):
            return new, transfers, extents, substratum, 0
        for rounds in range(1, 51):
            # A voxel is short only if more leaves than it holds and receives, so out > 0.
            short = new < 0
            available = y + into
            share = np.where(short, available * _MARGIN / np.where(short, out, 1.0), 1.0)
            # Below the smallest normal number no relative margin survives rounding:
            # scaled by 4/6, a transfer of one unit in the last place rounds back up
            # to one. A voxel holding that little stops giving instead.
            share = np.where(short & (available < _TINY), 0.0, share)
            transfers, extents, substratum = self._scale(
                transfers, extents, np.minimum(share, 1.0), substratum
            )
            new, into, out = self._apply(y, transfers, extents, substratum)
            if np.all(new >= 0):
                return new, transfers, extents, substratum, rounds
        over = out > y * _MARGIN  # the stock-only limiter
        share = np.where(over, y * _MARGIN / np.where(over, out, 1.0), 1.0)
        while True:
            transfers, extents, substratum = self._scale(transfers, extents, share, substratum)
            new, _, _ = self._apply(y, transfers, extents, substratum)
            if np.all(new >= 0):
                return new, transfers, extents, substratum, 51
            # Whatever rounding leaves below zero stops giving: a voxel that gives
            # nothing can only gain, so each pass settles more voxels for good.
            share = np.where(new < 0, 0.0, 1.0)

    # -- one step and an adaptive span ----------------------------------------------------

    def step(
        self, y: Field, h: float, reference: Field, atol: Tolerance, rtol: float
    ) -> tuple[Field, NDArray[np.float64], Field, int, int]:
        """One SDIRK2 step: new state, imports, error estimate, limiter rounds, Newton steps.

        Where the plaque spreads, half the step's wear comes off before it and
        half after, with the packing; what that detaches is subtracted from
        the imports, and kept in :attr:`_detaching` until the step is accepted.
        """
        # Traces far below any tolerance may underflow to zero; that loses nothing.
        with np.errstate(under="ignore"):
            tolerance = _tolerance(atol, y.ndim - 1)
            if self.spreading is None:
                return self._step(y, h, reference, tolerance, rtol)
            detached = np.zeros(y.shape[0])
            if self.spreading.wear_um_per_h > 0:
                y, before = self._detach(y, 0.5 * h)
                detached += before
            new, imports, estimate, rounds, iterations = self._step(
                y, h, reference, tolerance, rtol
            )
            new, after = self._detach(new, 0.5 * h)
            self._detaching = detached + after
            return new, imports - self._detaching, estimate, rounds, iterations

    def _detach(self, y: Field, hours: float) -> tuple[Field, NDArray[np.float64]]:
        """Pack the plaque, wear ``hours`` of it away; what passes its top leaves the box."""
        assert self.spreading is not None
        return self.spreading.project(y, hours)

    def _step(
        self, y: Field, h: float, reference: Field, atol: Tolerance, rtol: float
    ) -> tuple[Field, NDArray[np.float64], Field, int, int]:
        g = GAMMA
        reference = _magnitude(reference)
        system = self._system(y, g * h)
        stage1, system, n1 = self._stage(y, g * h, y, system, reference, atol, rtol)
        s1 = self.evaluate(stage1)
        base = y + (1 - g) * h * s1.rate
        stage2, system, n2 = self._stage(base, g * h, stage1, system, reference, atol, rtol)
        s2 = self.evaluate(stage2)
        transfers = tuple(
            h * ((1 - g) * f1 + g * f2) for f1, f2 in zip(s1.fluxes, s2.fluxes, strict=True)
        )
        extents = h * ((1 - g) * s1.reactions + g * s2.reactions)
        substratum = None
        if s1.substratum is not None and s2.substratum is not None:
            substratum = h * ((1 - g) * s1.substratum + g * s2.substratum)
        new, transfers, extents, substratum, rounds = self._limited_update(
            y, transfers, extents, substratum
        )
        scale = atol + rtol * np.maximum(reference, _magnitude(y, new))
        estimate, _ = system.solve(g * h * (s2.rate - s1.rate), scale=scale, tolerance=1e-2)
        top = -transfers[-1][..., -1]  # into the box through each top face
        imports = top.reshape(top.shape[0], -1).sum(axis=1) / self.spacing[-1]
        if substratum is not None:  # and through the substratum: bound less detached
            imports = imports + substratum.reshape(top.shape[0], -1).sum(axis=1) / self.spacing[-1]
        if self.air is not None:  # and from the air
            self._airing = self.air.exchanged(extents[self._air_rows])
            imports = imports + self._airing
        return new, imports, estimate, rounds, n1 + n2

    def starting_step(self, y: Field, reference: Field, atol: Tolerance, rtol: float) -> float:
        """A first step from the rates and their change over a trial Euler step.

        The algorithm of Hairer, Norsett and Wanner (1993, section II.4) for a
        second-order method, in the error norm of :meth:`integrate`: the Euler
        increment stays a hundredth of the state, and the local error, from the
        change in the rates, a hundredth of the tolerance.
        """
        with np.errstate(under="ignore"):
            return self._starting_step(y, reference, _tolerance(atol, y.ndim - 1), rtol)

    def _starting_step(self, y: Field, reference: Field, atol: Tolerance, rtol: float) -> float:
        scale = atol + rtol * np.maximum(_magnitude(reference), _magnitude(y))
        rate = self.evaluate(y).rate
        size = float(np.max(np.abs(y) / scale))
        speed = float(np.max(np.abs(rate) / scale))
        trial = 1e-6 if size < 1e-5 or speed < 1e-5 else 0.01 * size / speed
        change = self.evaluate(y + trial * rate).rate - rate
        curvature = float(np.max(np.abs(change) / scale)) / trial
        if max(speed, curvature) <= 1e-15:
            return max(1e-6, trial * 1e-3)
        return min(100 * trial, (0.01 / max(speed, curvature)) ** (1 / 3))

    def integrate(
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
        """Advance ``y`` by ``span`` in adaptive steps; imports are summed over the span.

        Without ``first_step``, the first step comes from :meth:`starting_step`.
        Imports are in the units the ledger sums: concentration summed over
        voxels, per component. ``absolute_tolerance`` is one number, or one per
        component. ``start_h`` is the time the span starts at, for a subclass
        whose rates depend on time (:mod:`marse.core.reservoir`); each step
        sees its own start as ``self._now``.
        """
        self._now = start_h
        absolute_tolerance = _tolerance(absolute_tolerance, y.ndim - 1)
        reference = _magnitude(y if peak is None else peak)
        if first_step is None:
            first_step = self.starting_step(y, reference, absolute_tolerance, relative_tolerance)
        if not span > 0 or not first_step > 0:
            raise ValueError("span and first_step must be positive")
        elapsed, h = 0.0, first_step
        imports = np.zeros(y.shape[0])
        accepted = rejected = limited = failures = newton = 0
        while span - elapsed > 1e-12 * span:
            step = min(h, span - elapsed)
            clipped = step < h
            self._now = start_h + elapsed
            try:
                new, entered, estimate, rounds, iterations = self.step(
                    y, step, reference, absolute_tolerance, relative_tolerance
                )
            except NewtonFailure:
                failures += 1
                rejected += 1
                h = step * 0.25
                if h <= 1e-14 * span:
                    raise
                continue
            newton += iterations
            scale = absolute_tolerance + relative_tolerance * np.maximum(
                reference, _magnitude(y, new)
            )
            with np.errstate(under="ignore"):
                error = float(np.max(np.abs(estimate) / scale))
            grow = 0.9 / math.sqrt(error) if error > 0 else 2.0
            if error <= 1.0:
                y = new
                imports += entered
                if self.spreading is not None:
                    self.detached += self._detaching
                if self.air is not None:
                    self.aired += self._airing
                elapsed += step
                accepted += 1
                limited += rounds > 0
                if not clipped:
                    h = step * min(2.0, max(0.2, grow))
            else:
                rejected += 1
                h = step * min(0.5, max(0.1, grow))
                if h <= 1e-14 * span:
                    raise ArithmeticError(
                        "the adaptive step shrank below 1e-14 of the span; the system is too "
                        "stiff for these tolerances"
                    )
        return y, imports, StepStats(accepted, rejected, limited, failures, newton, h)
