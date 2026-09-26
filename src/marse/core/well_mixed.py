"""The well-mixed engine: a reaction network running in a closed box.

``marse run`` sends a runnable version 2 configuration here (docs/networks.md).
Concentrations follow

    dc/dt = N^T r(c),

one rate per process driving every component that process touches
(docs/theory.md, sections 3.6 and 3.7). Consumption therefore follows growth
exactly, every species is updated together, and no process can create
matter: known defects 2 and 7 cannot arise here.

- **Integration** (:mod:`marse.core.integrators`) is positive and conservative
  by construction, in adaptive substeps under the configured tolerances.
- **The ledger** (:mod:`marse.core.ledger`) checks the carbon, nitrogen and
  electron totals after every step. A drift stops the run instead of
  returning it.
- **The manifest** records the configuration, the engine and integrator, the
  balance, the substeps taken, and a SHA-256 of the final state, so
  ``marse replay`` reproduces the run bit for bit.
"""

from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from marse.core.integrators import integrate
from marse.core.ledger import Ledger
from marse.core.provenance import Manifest
from marse.microbes.kinetics import compile_rates, process_rates
from marse.schemas.experiment import WellMixedConfig
from marse.schemas.formula import QUANTITIES

__all__ = [
    "ENGINE_VERSION",
    "INTEGRATOR_VERSION",
    "ExhaustedComponentError",
    "WellMixedResult",
    "run",
]

ENGINE_VERSION = "well_mixed_v2"
INTEGRATOR_VERSION = "limited_heun_adaptive_v1"


class ExhaustedComponentError(RuntimeError):
    """A component a process assumes to be in excess has run out.

    The configuration stated the assumption (``assumed_in_excess``). A run in
    which it fails is stopped rather than continued on a false premise.
    """


@dataclass(frozen=True, slots=True)
class WellMixedResult:
    """A finished run: the recorded trajectory, its balance and its manifest."""

    config: WellMixedConfig
    manifest: Manifest
    times_h: NDArray[np.float64]
    concentrations_mol_per_m3: NDArray[np.float64]  # shape (records, components)

    @property
    def component_names(self) -> tuple[str, ...]:
        return self.config.network.component_names

    @property
    def final_mol_per_m3(self) -> dict[str, float]:
        return dict(zip(self.component_names, self.concentrations_mol_per_m3[-1], strict=True))

    def write_trajectory(self, path: str | Path) -> Path:
        """Write the recorded trajectory as CSV, with units in the headers."""
        destination = Path(path)
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["time_h"] + [f"{n}_mol_per_m3" for n in self.component_names])
            for t, row in zip(self.times_h, self.concentrations_mol_per_m3, strict=True):
                writer.writerow([f"{t:.6f}"] + [f"{v:.10g}" for v in row])
        return destination


def _digest(names: tuple[str, ...], state: NDArray[np.float64]) -> str:
    """SHA-256 over the final state in a fixed byte order, independent of the machine."""
    digest = hashlib.sha256()
    for name, value in zip(names, np.asarray(state, dtype="<f8"), strict=True):
        digest.update(f"{name}:".encode())
        digest.update(value.tobytes())
    return digest.hexdigest()


def run(config: WellMixedConfig) -> WellMixedResult:
    """Integrate the network from its initial state to ``duration_h``.

    Raises :class:`~marse.core.simulation.ConservationError` if carbon,
    nitrogen or electrons drift, and :class:`ExhaustedComponentError` if a
    component a process assumes to be in excess runs out.
    """
    started = datetime.now(UTC)
    network = config.network
    names = network.component_names
    stoichiometry = network.stoichiometric_matrix()  # processes x components
    terms = compile_rates(network)
    index = {name: i for i, name in enumerate(names)}
    assumptions = [
        (p.name, [index[n] for n in p.rate.assumed_in_excess])
        for p in network.processes
        if p.rate is not None and p.rate.assumed_in_excess
    ]

    def rates(concentrations: NDArray[np.float64]) -> NDArray[np.float64]:
        return process_rates(terms, concentrations)

    state = np.array([config.initial_mol_per_m3[n] for n in names], dtype=float)
    ledger = Ledger(network.composition_matrix(), state, QUANTITIES)
    peak = state.copy()
    times, rows = [0.0], [state.copy()]
    now, substep = 0.0, config.timestep_h
    accepted = rejected = limited = 0
    steps = config.steps

    def assumptions_hold(concentrations: NDArray[np.float64]) -> None:
        for process, indices in assumptions:
            for i in indices:
                if concentrations[i] <= config.absolute_tolerance_mol_per_m3:
                    raise ExhaustedComponentError(
                        f"'{names[i]}' ran out in the step from t = {now:.6g} h, but process "
                        f"'{process}' assumes it is in excess; give that process a monod "
                        f"factor for '{names[i]}' instead"
                    )

    for step in range(1, steps + 1):
        if assumptions:
            assumptions_hold(state)
        # Times come from the step count, so they never accumulate rounding.
        target = min(step * config.timestep_h, config.duration_h)
        state, stats = integrate(
            state,
            target - now,
            rates,
            stoichiometry,
            first_step=substep,
            relative_tolerance=config.relative_tolerance,
            absolute_tolerance=config.absolute_tolerance_mol_per_m3,
            peak=peak,
            guard=assumptions_hold if assumptions else None,
        )
        substep, now = stats.next_step, target
        accepted += stats.accepted
        rejected += stats.rejected
        limited += stats.limited
        peak = np.maximum(peak, state)
        ledger.check(state, step=step, time_h=now)
        if step % config.record_every == 0 or step == steps:
            times.append(now)
            rows.append(state.copy())

    outputs: dict[str, Any] = {
        "final_time_h": now,
        "final_mol_per_m3": {n: float(v) for n, v in zip(names, state, strict=True)},
        "balance": ledger.summary(state),
        "substeps": {"accepted": accepted, "rejected": rejected, "limited": limited},
        "final_state_sha256": _digest(names, state),
    }
    manifest = Manifest.build(
        config=config,
        models={"engine": ENGINE_VERSION, "integrator": INTEGRATOR_VERSION},
        random_streams=(),
        started=started,
        finished=datetime.now(UTC),
        steps=steps,
        outputs=outputs,
    )
    return WellMixedResult(
        config=config,
        manifest=manifest,
        times_h=np.array(times),
        concentrations_mol_per_m3=np.array(rows),
    )
