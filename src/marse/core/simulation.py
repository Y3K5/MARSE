"""The deterministic simulation kernel.

Phase 1 of the roadmap: a well-mixed batch culture of one or more organisms
competing for a single limiting substrate. There is no space here — that is
Phase 2. What this module owns is everything that is *not* biology: the clock,
the canonical state, the order in which models are applied, seeded randomness,
checkpoints, conservation checks and the run manifest.

The biology comes from the tested providers in :mod:`marse.microbes`, which
this module calls without reimplementing:

    dX_i/dt = (mu_i - b_i) X_i
    dS/dt   = -sum_i (mu_i / Y_i + m_i) X_i
    mu_i    = mu_opt,i * gamma_T(T) * gamma_pH(pH) * S / (K_S,i + S)

Integration uses classical fourth-order Runge-Kutta at a fixed step, so a run
is reproducible to the bit from its configuration (docs/theory.md, section 9).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from marse.core.config import ExperimentConfig
from marse.core.provenance import Manifest
from marse.core.seeds import SeedRegistry
from marse.core.state import SimulationState
from marse.microbes.cardinal import cardinal_ph, cardinal_temperature
from marse.microbes.growth import monod

__all__ = ["ConservationError", "SimulationResult", "UnstableStepError", "run"]

MODELS = {
    "growth": "monod_v1",
    "uptake": "pirt_v1",
    "temperature": "ctmi_v1",
    "ph": "cpm_v1",
    "integrator": "rk4_v1",
}
"""Versioned provider names, recorded in the manifest so a result names its models."""

CONSERVATION_TOLERANCE_MM = 1e-9
"""Absolute tolerance on the substrate balance, checked every step."""

_NEGATIVE_BIOMASS_TOLERANCE = 1e-9
"""Relative size of a negative biomass excursion that counts as instability."""


class ConservationError(RuntimeError):
    """Substrate was created or destroyed, or the state stopped being finite."""


class UnstableStepError(RuntimeError):
    """The timestep is too large for the configured rates."""


@dataclass(frozen=True, slots=True)
class SimulationResult:
    """Outcome of a run: the trajectory, the final state, and the manifest."""

    config: ExperimentConfig
    manifest: Manifest
    times_h: NDArray[np.float64]
    biomass_g_per_l: NDArray[np.float64]  # shape (checkpoints, organisms)
    substrate_mm: NDArray[np.float64]
    final_state: SimulationState

    @property
    def organism_names(self) -> tuple[str, ...]:
        return tuple(o.name for o in self.config.organisms)

    def write_trajectory(self, path: str | Path) -> Path:
        """Write the checkpointed trajectory as CSV, with units in the headers."""
        destination = Path(path)
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                ["time_h", f"{self.config.substrate.name}_mm"]
                + [f"{name}_g_per_l" for name in self.organism_names]
            )
            for i, t in enumerate(self.times_h):
                writer.writerow(
                    [f"{t:.6f}", f"{self.substrate_mm[i]:.10g}"]
                    + [f"{value:.10g}" for value in self.biomass_g_per_l[i]]
                )
        return destination


def _max_growth_rates(config: ExperimentConfig) -> NDArray[np.float64]:
    """Scale each organism's optimal rate by the environment (theory.md section 2)."""
    temperature = config.environment.temperature_c
    ph = config.environment.ph
    rates = np.empty(len(config.organisms), dtype=float)
    for i, organism in enumerate(config.organisms):
        gamma = 1.0
        if organism.cardinal_temperature_c is not None:
            gamma *= float(cardinal_temperature(temperature, *organism.cardinal_temperature_c))
        if organism.cardinal_ph is not None:
            gamma *= float(cardinal_ph(ph, *organism.cardinal_ph))
        rates[i] = organism.mu_opt_per_h * gamma
    return rates


def _derivatives(
    biomass: NDArray[np.float64],
    substrate: float,
    mu_max: NDArray[np.float64],
    k_s: NDArray[np.float64],
    inverse_yield: NDArray[np.float64],
    maintenance: NDArray[np.float64],
    decay: NDArray[np.float64],
) -> tuple[NDArray[np.float64], float]:
    """Return (dX/dt, consumption rate). Substrate follows from the consumption."""
    live = np.maximum(biomass, 0.0)
    mu = np.asarray(monod(max(substrate, 0.0), mu_max, k_s), dtype=float)
    growth = (mu - decay) * live
    consumption = float(np.sum((mu * inverse_yield + maintenance) * live))
    return growth, consumption


def run(config: ExperimentConfig) -> SimulationResult:
    """Run a batch simulation to completion.

    Raises :class:`ConservationError` if the substrate balance drifts, rather
    than returning a result that silently created or destroyed mass.
    """
    started = datetime.now(UTC)
    seeds = SeedRegistry(config.seed)
    # Phase 1 is deterministic; the streams exist so that stochastic providers
    # added in Phase 3 inherit reproducibility rather than bolting it on.
    seeds.stream("core")

    names = tuple(o.name for o in config.organisms)
    mu_max = _max_growth_rates(config)
    k_s = np.array([o.k_s_mm for o in config.organisms])
    inverse_yield = np.array([1.0 / o.yield_g_per_mmol for o in config.organisms])
    maintenance = np.array([o.maintenance_mmol_per_g_per_h for o in config.organisms])
    decay = np.array([o.decay_per_h for o in config.organisms])

    state = SimulationState.initial(
        biomass=np.array([o.initial_biomass_g_per_l for o in config.organisms]),
        substrate_mm=config.substrate.initial_mm,
    )
    initial_substrate = state.substrate_mm

    times = [state.time_h]
    biomass_history = [state.biomass_g_per_l.copy()]
    substrate_history = [state.substrate_mm]
    next_checkpoint = config.checkpoint_interval_h

    total_steps = config.steps
    for step in range(1, total_steps + 1):
        # The final step is shortened so the run lands exactly on duration_h.
        t0 = state.time_h
        dt = min(config.timestep_h, config.duration_h - t0)
        if dt <= 0.0:
            break

        def rates(
            biomass: NDArray[np.float64], substrate: float
        ) -> tuple[NDArray[np.float64], float]:
            return _derivatives(biomass, substrate, mu_max, k_s, inverse_yield, maintenance, decay)

        x0, s0 = state.biomass_g_per_l, state.substrate_mm
        # A diverging run produces overflow and invalid operations. Those are
        # detected explicitly below, which is deterministic; relying on NumPy's
        # warnings would make the behaviour depend on the caller's warning filters.
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            k1x, k1c = rates(x0, s0)
            k2x, k2c = rates(x0 + dt / 2 * k1x, s0 - dt / 2 * k1c)
            k3x, k3c = rates(x0 + dt / 2 * k2x, s0 - dt / 2 * k2c)
            k4x, k4c = rates(x0 + dt * k3x, s0 - dt * k3c)

            biomass = x0 + dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
            consumed = dt / 6 * (k1c + 2 * k2c + 2 * k3c + k4c)

        if not (np.all(np.isfinite(biomass)) and np.isfinite(consumed)):
            raise ConservationError(
                f"the state became non-finite at step {step} (t = {t0 + dt:.4f} h); "
                "the run diverged"
            )

        # A fixed-step explicit method can overshoot if the step is too large for
        # the dynamics. Biomass driven appreciably below zero is the signature, and
        # clamping it silently would turn an unstable run into a plausible-looking
        # wrong answer, so it is reported instead.
        most_negative = float(np.min(biomass)) if biomass.size else 0.0
        if most_negative < -_NEGATIVE_BIOMASS_TOLERANCE * max(1.0, float(np.max(np.abs(x0)))):
            raise UnstableStepError(
                f"biomass reached {most_negative:.3e} g/L at step {step} (t = {t0 + dt:.4f} h): "
                f"timestep_h = {config.timestep_h:g} is too large for these rates; reduce it"
            )

        # Substrate cannot go below zero; clamp and consume only what was there.
        consumed = min(consumed, state.substrate_mm)
        substrate = state.substrate_mm - consumed
        biomass = np.maximum(biomass, 0.0)

        state = state.advanced(
            time_h=t0 + dt,
            step=step,
            biomass=biomass,
            substrate_mm=substrate,
            substrate_consumed_mm=state.substrate_consumed_mm + consumed,
        )

        # Runtime invariant: every millimole is either still dissolved or consumed.
        balance = state.substrate_mm + state.substrate_consumed_mm - initial_substrate
        if not np.isfinite(balance) or abs(balance) > CONSERVATION_TOLERANCE_MM:
            raise ConservationError(
                f"substrate balance drifted by {balance:.3e} mM at step {step} "
                f"(t = {state.time_h:.4f} h); the run is not mass-conserving"
            )

        if state.time_h >= next_checkpoint - 1e-12 or step == total_steps:
            times.append(state.time_h)
            biomass_history.append(state.biomass_g_per_l.copy())
            substrate_history.append(state.substrate_mm)
            while next_checkpoint <= state.time_h + 1e-12:
                next_checkpoint += config.checkpoint_interval_h

    finished = datetime.now(UTC)
    final = state.to_dict(names)
    manifest = Manifest.build(
        config=config,
        models=MODELS,
        random_streams=seeds.issued_names,
        started=started,
        finished=finished,
        steps=state.step,
        outputs={
            "final_state": final,
            "substrate_consumed_mm": state.substrate_consumed_mm,
            "checkpoints": len(times),
            "max_growth_rate_per_h": dict(zip(names, (float(v) for v in mu_max), strict=True)),
        },
    )
    return SimulationResult(
        config=config,
        manifest=manifest,
        times_h=np.array(times),
        biomass_g_per_l=np.array(biomass_history),
        substrate_mm=np.array(substrate_history),
        final_state=state,
    )
