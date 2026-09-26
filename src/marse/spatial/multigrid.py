"""Geometric multigrid for the linear systems of implicit reaction-diffusion steps.

An implicit step (docs/theory.md, section 9.7) needs, at every Newton
iteration, the solution of

    x - a (L x + B x) = b,

for ``x`` of shape ``(components, *shape)``, where

- ``L`` is diffusion (:mod:`marse.spatial.transport`) with the bulk held at
  zero, as a correction equation requires;
- ``B`` holds one small matrix per voxel, the Jacobian of the reactions there,
  which couples the components within a voxel but not across voxels;
- ``a`` is the step weight, gamma h.

Diffusion couples voxels over the whole grid, and with the long steps an
implicit scheme exists to take, ``a D / h^2`` reaches 10^4 or more. Simple
iterations then stall on the smooth part of the error. Multigrid
(Briggs, Henson and McCormick 2000) removes that part on coarser grids, where
it is no longer smooth, so each cycle reduces the error by a similar factor on
any grid. Its ingredients here:

- **smoothing**: red-black Gauss-Seidel over voxels, each voxel's components
  solved together through its small block;
- **restriction**: averaging the children of a coarse voxel;
- **prolongation**: linear interpolation between voxel centres;
- **coarse operators**: re-discretised on each level, with the reaction blocks
  averaged over children;
- **the coarsest level**: solved exactly.

The cycle preconditions a restarted GMRES (Saad and Schultz 1986), which keeps
convergence robust when the reaction blocks make the system non-symmetric.
Everything is deterministic: the same system gives the same answer, bit for
bit.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse.spatial.transport import diagonal, divergence, face_fluxes

__all__ = ["ImplicitSystem", "gmres"]

_DIRECT_LIMIT = 512  # coarsest levels up to this many unknowns are solved exactly


@dataclass(slots=True)
class _Level:
    shape: tuple[int, ...]
    spacing: tuple[float, ...]
    blocks: NDArray[np.float64]  # (voxels, J, J)
    diag_l: NDArray[np.float64]  # (J, voxels)
    inverse: NDArray[np.float64]  # (voxels, J, J), of the diagonal blocks of the system
    colours: tuple[tuple[NDArray[np.intp], NDArray[np.float64]], ...]  # voxels, their inverses
    couplings: tuple[
        NDArray[np.float64], ...
    ]  # per axis, D / h^2 per component, shaped to broadcast
    coarsened: tuple[int, ...]  # axes halved to reach the next level
    dense_inverse: NDArray[np.float64] | None = None


class ImplicitSystem:
    """The system x - a (L x + B x) = b on a grid, with its multigrid hierarchy.

    ``blocks`` has shape ``(J, J, *shape)``: at each voxel, the derivative of the
    reaction rates of every component with respect to every component.
    """

    def __init__(
        self,
        shape: tuple[int, ...],
        spacing: tuple[float, ...],
        diffusivity: NDArray[np.float64],
        a: float,
        blocks: NDArray[np.float64],
    ) -> None:
        self.components = int(diffusivity.size)
        self.diffusivity = np.asarray(diffusivity, dtype=float)
        self.a = float(a)
        voxel_blocks = np.moveaxis(blocks.reshape(self.components, self.components, -1), -1, 0)
        self.levels: list[_Level] = []
        level = self._make_level(tuple(shape), tuple(spacing), voxel_blocks)
        while True:
            self.levels.append(level)
            axes = tuple(axis for axis, n in enumerate(level.shape) if n % 2 == 0 and n >= 2)
            voxels = math.prod(level.shape)
            if not axes or voxels * self.components <= _DIRECT_LIMIT:
                break
            level.coarsened = axes
            shape_c = tuple(n // 2 if axis in axes else n for axis, n in enumerate(level.shape))
            spacing_c = tuple(h * 2 if axis in axes else h for axis, h in enumerate(level.spacing))
            blocks_c = _restrict_blocks(level.blocks, level.shape, axes)
            level = self._make_level(shape_c, spacing_c, blocks_c)
        coarsest = self.levels[-1]
        unknowns = math.prod(coarsest.shape) * self.components
        if unknowns <= _DIRECT_LIMIT:
            coarsest.dense_inverse = np.linalg.inv(self._dense(coarsest))

    # -- construction -----------------------------------------------------------------

    def _make_level(
        self, shape: tuple[int, ...], spacing: tuple[float, ...], blocks: NDArray[np.float64]
    ) -> _Level:
        j = self.components
        diag_l = diagonal(self.diffusivity, spacing, shape).reshape(j, -1)
        system = np.eye(j)[None] - self.a * blocks
        system[:, np.arange(j), np.arange(j)] -= self.a * diag_l.T
        inverse = np.linalg.inv(system)
        parity = (np.indices(shape).sum(axis=0) % 2 == 0).reshape(-1)
        colours = tuple(
            (index, np.ascontiguousarray(inverse[index]))
            for index in (np.flatnonzero(parity), np.flatnonzero(~parity))
        )
        dims = len(shape)
        couplings = tuple((self.diffusivity / h**2).reshape((-1,) + (1,) * dims) for h in spacing)
        return _Level(
            shape=shape,
            spacing=spacing,
            blocks=blocks,
            diag_l=diag_l,
            inverse=inverse,
            colours=colours,
            couplings=couplings,
            coarsened=(),
        )

    def _dense(self, level: _Level) -> NDArray[np.float64]:
        n = math.prod(level.shape) * self.components
        matrix = np.empty((n, n))
        unit = np.zeros(n)
        for k in range(n):
            unit[k] = 1.0
            matrix[:, k] = self._apply(level, unit.reshape(self.components, -1)).reshape(-1)
            unit[k] = 0.0
        return matrix

    # -- operators on one level ---------------------------------------------------------

    def _diffuse(self, level: _Level, x: NDArray[np.float64]) -> NDArray[np.float64]:
        field = x.reshape(self.components, *level.shape)
        rate = divergence(face_fluxes(field, self.diffusivity, level.spacing, None), level.spacing)
        return rate.reshape(self.components, -1)

    def _apply(self, level: _Level, x: NDArray[np.float64]) -> NDArray[np.float64]:
        react = np.einsum("vjk,kv->jv", level.blocks, x)
        return x - self.a * (self._diffuse(level, x) + react)

    def _neighbours(self, level: _Level, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """The off-diagonal part of diffusion: each voxel's neighbours, weighted by D / h^2."""
        field = x.reshape(self.components, *level.shape)
        total = np.zeros_like(field)
        dims = len(level.shape)
        for axis in range(dims - 1):
            if level.shape[axis] >= 2:
                total += level.couplings[axis] * (
                    np.roll(field, 1, axis=axis + 1) + np.roll(field, -1, axis=axis + 1)
                )
        coupling = level.couplings[-1]
        total[..., 1:] += coupling * field[..., :-1]
        total[..., :-1] += coupling * field[..., 1:]
        return total.reshape(self.components, -1)

    def _smooth(
        self, level: _Level, x: NDArray[np.float64], b: NDArray[np.float64], sweeps: int
    ) -> NDArray[np.float64]:
        for _ in range(sweeps):
            for index, inverse in level.colours:
                rhs = b[:, index] + self.a * self._neighbours(level, x)[:, index]
                x[:, index] = np.einsum("vjk,kv->jv", inverse, rhs)
        return x

    def _cycle(self, index: int, b: NDArray[np.float64]) -> NDArray[np.float64]:
        level = self.levels[index]
        if index == len(self.levels) - 1:
            if level.dense_inverse is not None:
                return (level.dense_inverse @ b.reshape(-1)).reshape(b.shape)
            return self._smooth(level, np.zeros_like(b), b, 40)
        x = self._smooth(level, np.zeros_like(b), b, 2)
        residual = b - self._apply(level, x)
        coarse = _restrict(residual, level.shape, level.coarsened, self.components)
        correction = self._cycle(index + 1, coarse)
        x += _prolong(
            correction,
            self.levels[index + 1].shape,
            level.coarsened,
            len(level.shape),
            self.components,
        )
        return self._smooth(level, x, b, 2)

    # -- public -------------------------------------------------------------------------

    def apply(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """The system applied to x, shape (components, *shape)."""
        return self._apply(self.levels[0], x.reshape(self.components, -1)).reshape(x.shape)

    def precondition(self, b: NDArray[np.float64]) -> NDArray[np.float64]:
        """One multigrid V-cycle from zero: an approximate solution of the system for b."""
        return self._cycle(0, b.reshape(self.components, -1)).reshape(b.shape)

    def solve(
        self,
        b: NDArray[np.float64],
        *,
        scale: NDArray[np.float64],
        tolerance: float = 1e-6,
        max_iterations: int = 200,
    ) -> tuple[NDArray[np.float64], int]:
        """Solve to ``tolerance`` in the norm weighted by 1/scale; return x and iterations."""
        return gmres(
            self.apply,
            b,
            self.precondition,
            scale=scale,
            tolerance=tolerance,
            max_iterations=max_iterations,
        )


def gmres(
    apply: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    b: NDArray[np.float64],
    precondition: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    *,
    scale: NDArray[np.float64],
    tolerance: float,
    restart: int = 20,
    max_iterations: int = 200,
) -> tuple[NDArray[np.float64], int]:
    """Right-preconditioned restarted GMRES in the norm weighted by 1/scale.

    Components differ in size by many orders of magnitude (biomass near 10^3
    mol/m3, oxygen near 10^-1), so residuals are measured relative to the same
    per-entry scale the error control uses.
    """
    weight = 1.0 / scale

    def norm(v: NDArray[np.float64]) -> float:
        return float(np.sqrt(np.sum((v * weight) ** 2)))

    def dot(u: NDArray[np.float64], v: NDArray[np.float64]) -> float:
        return float(np.sum(u * v * weight**2))

    x = np.zeros_like(b)
    target = tolerance * norm(b)
    residual = b.copy()
    beta = norm(residual)
    used = 0
    if beta <= target or beta == 0.0:
        return x, 0
    while used < max_iterations:
        basis = [residual / beta]
        directions = []
        hessenberg = np.zeros((restart + 1, restart))
        coefficients = np.zeros(0)
        for k in range(restart):
            z = precondition(basis[k])
            directions.append(z)
            w = apply(z)
            for i in range(k + 1):
                hessenberg[i, k] = dot(w, basis[i])
                w = w - hessenberg[i, k] * basis[i]
            hessenberg[k + 1, k] = norm(w)
            used += 1
            rhs = np.zeros(k + 2)
            rhs[0] = beta
            coefficients, *_ = np.linalg.lstsq(hessenberg[: k + 2, : k + 1], rhs, rcond=None)
            achieved = float(np.linalg.norm(rhs - hessenberg[: k + 2, : k + 1] @ coefficients))
            if achieved <= target or hessenberg[k + 1, k] == 0.0 or used >= max_iterations:
                break
            basis.append(w / hessenberg[k + 1, k])
        for c, z in zip(coefficients, directions, strict=False):
            x = x + c * z
        residual = b - apply(x)
        beta = norm(residual)
        if beta <= target:
            return x, used
    raise ArithmeticError(
        f"the linear solver did not converge in {max_iterations} iterations "
        f"(residual {beta / max(norm(b), 1e-300):.2e} of the right-hand side)"
    )


def _restrict(
    x: NDArray[np.float64], shape: tuple[int, ...], axes: tuple[int, ...], j: int
) -> NDArray[np.float64]:
    field = x.reshape(j, *shape)
    for axis in axes:
        a = axis + 1
        even = np.take(field, np.arange(0, field.shape[a], 2), axis=a)
        odd = np.take(field, np.arange(1, field.shape[a], 2), axis=a)
        field = 0.5 * (even + odd)
    return field.reshape(j, -1)


def _restrict_blocks(
    blocks: NDArray[np.float64], shape: tuple[int, ...], axes: tuple[int, ...]
) -> NDArray[np.float64]:
    j = blocks.shape[-1]
    field = blocks.reshape(*shape, j, j)
    for axis in axes:
        even = np.take(field, np.arange(0, field.shape[axis], 2), axis=axis)
        odd = np.take(field, np.arange(1, field.shape[axis], 2), axis=axis)
        field = 0.5 * (even + odd)
    return field.reshape(-1, j, j)


def _prolong(
    x: NDArray[np.float64],
    coarse_shape: tuple[int, ...],
    axes: tuple[int, ...],
    dims: int,
    j: int,
) -> NDArray[np.float64]:
    """Linear interpolation from coarse voxel centres to fine ones, with the boundary conditions."""
    field = x.reshape(j, *coarse_shape)
    for axis in axes:
        a = axis + 1
        if axis < dims - 1:  # periodic
            below = np.roll(field, 1, axis=a)
            above = np.roll(field, -1, axis=a)
        else:  # height: no flux at the substratum, zero correction at the top face
            first = np.take(field, [0], axis=a)
            last = np.take(field, [field.shape[a] - 1], axis=a)
            below = np.concatenate(
                (first, np.take(field, np.arange(field.shape[a] - 1), axis=a)), axis=a
            )
            above = np.concatenate(
                (np.take(field, np.arange(1, field.shape[a]), axis=a), -last), axis=a
            )
        even = 0.75 * field + 0.25 * below
        odd = 0.75 * field + 0.25 * above
        stacked = np.stack((even, odd), axis=a + 1)
        new_shape = list(field.shape)
        new_shape[a] *= 2
        field = stacked.reshape(new_shape)
    return field.reshape(j, -1)
