"""A direct solve for one-dimensional columns: block-tridiagonal elimination.

In a column, the implicit system of :mod:`marse.spatial.multigrid`,

    x - a (L x + B x) = b,

couples each voxel only to the one below and the one above. Diffusion
couples a component only to itself, so the blocks off the diagonal are
diagonal matrices, a D / h^2 per component; the J x J reaction block B sits
on the diagonal. Block Gaussian elimination from the substratum up, then
back substitution from the top (the Thomas algorithm, in blocks), solves the
system exactly in O(n J^3). Multigrid, which a column does not need, would
iterate towards the same answer.

The matrix is dominated by its diagonal wherever diffusion or a short step
makes it so, and it is the matrix of a stable implicit step elsewhere; the
eliminated blocks stay well conditioned in both cases (docs/theory.md,
section 9.7).

Each eliminated block is inverted by LAPACK, which factorises a matrix of
fewer than 100 unknowns on one thread, so a column replays bit for bit on any
number of threads. A network of 100 components or more uses
:func:`marse.spatial.multigrid.inverse`, which guarantees the same at any
size but is slower on small blocks.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from marse.spatial.multigrid import inverse
from marse.spatial.transport import diagonal, divergence, face_fluxes

__all__ = ["ColumnSystem"]

_SINGLE_THREADED = 100  # LAPACK factorises smaller matrices on one thread


class ColumnSystem:
    """The system x - a (L x + B x) = b on a column of voxels, solved directly.

    It has the interface of :class:`marse.spatial.multigrid.ImplicitSystem`.
    ``blocks`` has shape ``(J, J, n)``: at each voxel, the derivative of every
    reaction rate with respect to every component.
    """

    def __init__(
        self,
        shape: tuple[int, ...],
        spacing: tuple[float, ...],
        diffusivity: NDArray[np.float64],
        a: float,
        blocks: NDArray[np.float64],
        *,
        closed_top: bool = False,
    ) -> None:
        if len(shape) != 1:
            raise ValueError(f"a column has one axis, got the shape {shape}")
        (n,) = shape
        j = int(diffusivity.size)
        self.shape, self.spacing, self.components = tuple(shape), tuple(spacing), j
        self.diffusivity = np.asarray(diffusivity, dtype=float)
        self.a = float(a)
        self.closed_top = closed_top
        self.blocks = np.moveaxis(np.asarray(blocks, dtype=float).reshape(j, j, n), -1, 0)
        own = diagonal(self.diffusivity, self.spacing, self.shape, closed_top=closed_top)
        self.coupling = -self.a * self.diffusivity / self.spacing[-1] ** 2  # to each neighbour
        system = np.eye(j)[None] - self.a * self.blocks
        system[:, np.arange(j), np.arange(j)] -= self.a * own.reshape(j, n).T
        # Eliminate downwards: each voxel's block less what the one below passes up.
        pair = self.coupling[:, None] * self.coupling[None, :]
        invert = np.linalg.inv if j < _SINGLE_THREADED else inverse
        eliminated = np.empty((n, j, j))
        eliminated[0] = invert(system[0])
        for i in range(1, n):
            eliminated[i] = invert(system[i] - pair * eliminated[i - 1])
        self._eliminated = eliminated

    def apply(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """The system applied to x, shape (components, n)."""
        field = x.reshape(self.components, *self.shape)
        fluxes = face_fluxes(
            field, self.diffusivity, self.spacing, None, closed_top=self.closed_top
        )
        react = np.einsum("vjk,kv->jv", self.blocks, field)
        return (field - self.a * (divergence(fluxes, self.spacing) + react)).reshape(x.shape)

    def precondition(self, b: NDArray[np.float64]) -> NDArray[np.float64]:
        """The exact solution of the system for b."""
        rhs = b.reshape(self.components, -1).T  # (n, J)
        n = rhs.shape[0]
        forward = np.empty_like(rhs)
        forward[0] = self._eliminated[0] @ rhs[0]
        for i in range(1, n):
            forward[i] = self._eliminated[i] @ (rhs[i] - self.coupling * forward[i - 1])
        x = np.empty_like(rhs)
        x[-1] = forward[-1]
        for i in range(n - 2, -1, -1):
            x[i] = forward[i] - self._eliminated[i] @ (self.coupling * x[i + 1])
        return x.T.reshape(b.shape)

    def solve(
        self,
        b: NDArray[np.float64],
        *,
        scale: NDArray[np.float64],
        tolerance: float = 1e-6,
        max_iterations: int = 200,
    ) -> tuple[NDArray[np.float64], int]:
        """The exact solution, and one "iteration", for the interface of the multigrid solver."""
        del scale, tolerance, max_iterations  # exact: nothing to converge
        return self.precondition(b), 1
