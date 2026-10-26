"""The balance of carbon, nitrogen and electrons, checked as a run proceeds.

For each conserved quantity k, with composition I (components by quantities),

    residual_k(t) = M_k(t) - M_k(0) - imports_k + exports_k,  M_k = sum_j I_jk c_j,

must stay at rounding level (docs/theory.md, section 9.5). In a closed box
imports and exports are zero; in space, :meth:`Ledger.exchange` books what
crossed the top face, computed from the face fluxes independently of the
stored amounts. A residual
beyond the tolerance stops the run with a
:class:`~marse.core.simulation.ConservationError`, because a result that
created or destroyed matter must not be returned as if it were valid. The
largest residual seen is reported with every run, so each run proves its own
balance.

The tolerance is relative to the total absolute content, sum_j |I_jk| c_j(0),
which is never zero for a quantity the run contains. It sits far above
rounding, about 1e-15 per step, and far below any real leak.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from marse.core.simulation import ConservationError

__all__ = ["Ledger"]

DEFAULT_TOLERANCE = 1e-9


class Ledger:
    """Track conserved totals from an initial state and refuse any drift beyond tolerance."""

    def __init__(
        self,
        composition: NDArray[np.float64],
        initial: NDArray[np.float64],
        quantities: tuple[str, ...],
        tolerance: float = DEFAULT_TOLERANCE,
    ) -> None:
        self.composition = np.asarray(composition, dtype=float)
        self.quantities = quantities
        self.tolerance = tolerance
        self.initial = self.totals(initial)
        content = np.abs(self.composition).T @ np.asarray(initial, dtype=float).reshape(
            len(self.composition), -1
        ).sum(axis=1)
        self._content = content
        self._scale = np.where(content > 0, content, 1.0)
        self._exchanged = np.zeros(len(quantities))
        self._open = False
        self.imported = np.zeros(len(self.composition))
        self.largest_relative_residual = np.zeros(len(quantities))

    def exchange(self, imports: NDArray[np.float64]) -> None:
        """Book what crossed the boundary, per component: positive in, negative out.

        Imports are in the units of the state summed over cells. The tolerance
        grows with everything that has crossed, so a long open run is held to
        the same relative standard as a closed one.
        """
        amounts = np.asarray(imports, dtype=float)
        self._open = True
        self.imported = self.imported + amounts
        self._exchanged = self._exchanged + np.abs(self.composition).T @ np.abs(amounts)
        content = self._content + self._exchanged
        self._scale = np.where(content > 0, content, 1.0)

    def totals(self, state: NDArray[np.float64]) -> NDArray[np.float64]:
        """Each quantity summed over components (and cells, if any)."""
        flat = np.asarray(state, dtype=float).reshape(len(self.composition), -1)
        return self.composition.T @ flat.sum(axis=1)

    def check(self, state: NDArray[np.float64], *, step: int, time_h: float) -> None:
        """Record the residual now, and raise if any quantity has drifted too far."""
        residual = self.totals(state) - self.initial - self.composition.T @ self.imported
        relative = np.abs(residual) / self._scale
        self.largest_relative_residual = np.maximum(self.largest_relative_residual, relative)
        worst = int(np.argmax(relative))
        if not relative[worst] <= self.tolerance:  # also catches NaN
            raise ConservationError(
                f"{self.quantities[worst]} is not conserved at step {step} (t = {time_h:.6g} h): "
                f"its total changed by {residual[worst]:.3e}, {relative[worst]:.1e} of the "
                f"total content, beyond the tolerance of {self.tolerance:g}"
            )

    def summary(self, state: NDArray[np.float64]) -> dict[str, dict[str, float]]:
        """Per quantity: initial and final totals, imports if open, and the largest residual."""
        final = self.totals(state)
        imported = self.composition.T @ self.imported
        summary: dict[str, dict[str, float]] = {}
        for k, quantity in enumerate(self.quantities):
            entry = {"initial": float(self.initial[k]), "final": float(final[k])}
            if self._open:
                entry["imported"] = float(imported[k])
            entry["largest_relative_residual"] = float(self.largest_relative_residual[k])
            summary[quantity] = entry
        return summary
