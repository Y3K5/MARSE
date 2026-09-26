"""The spatial engine: a reaction network in a box of voxels over a surface.

``marse run`` sends a version 2 configuration with a ``domain`` here
(docs/networks.md, "Running a network in space"). Concentrations change as

    dc/dt = diffusion + N^T r(c),

in every voxel at once, integrated implicitly and conservatively by
:mod:`marse.core.implicit` (docs/theory.md, sections 4.7 and 9.8). Solutes
diffuse at physical diffusivities; biomass grows where it is, until Stage 2d
lets colonies spread.

- **The ledger** (:mod:`marse.core.ledger`) books what crosses the top face,
  computed from the face fluxes, and checks after every step that the carbon,
  nitrogen and electrons in the box changed by exactly that.
- **Frames** of every component's field go, as they are recorded, to a caller's
  function, which writes them to disk (``marse run`` writes a frame store and
  VTK files for ParaView).
- **The manifest** records the configuration, the engine and integrator, the
  balance, the steps taken, the time scales, and a SHA-256 of the final fields,
  so ``marse replay`` reproduces the run bit for bit.

Amounts are reported per m² of substratum, the unit biofilm measurements use:
a 1-D column and a 3-D box of the same chemistry report the same numbers.
"""

from __future__ import annotations

import csv
import hashlib
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from marse.core.implicit import ReactionTransport
from marse.core.ledger import Ledger
from marse.core.provenance import Manifest
from marse.microbes.kinetics import compile_rates, process_rates, rate_jacobian
from marse.schemas.experiment import ReactiveTransportConfig
from marse.schemas.formula import QUANTITIES
from marse.spatial.transport import Diffusion

__all__ = ["ENGINE_VERSION", "INTEGRATOR_VERSION", "ReactiveTransportResult", "run"]

ENGINE_VERSION = "reactive_transport_v2"
INTEGRATOR_VERSION = "sdirk2_limited_v1"

type FrameWriter = Callable[[int, float, NDArray[np.float64]], None]

_AMOL_PER_UM2_TO_MOL_PER_M2 = 1e-6  # 1e-18 mol per 1e-12 m2


@dataclass(frozen=True, slots=True)
class ReactiveTransportResult:
    """A finished spatial run: areal totals over time, the final fields, the manifest."""

    config: ReactiveTransportConfig
    manifest: Manifest
    times_h: NDArray[np.float64]
    totals_mol_per_m2: NDArray[np.float64]  # (records, components)
    imported_mol_per_m2: NDArray[np.float64]  # (records, components), since the start
    final_state: NDArray[np.float64]  # (components, *voxels)

    @property
    def component_names(self) -> tuple[str, ...]:
        return self.config.network.component_names

    def write_totals(self, path: str | Path) -> Path:
        """Write the areal totals and what has crossed the top face, as CSV with units."""
        destination = Path(path)
        names = self.component_names
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                ["time_h"]
                + [f"{n}_mol_per_m2" for n in names]
                + [f"{n}_imported_mol_per_m2" for n in names]
            )
            for t, totals, imported in zip(
                self.times_h, self.totals_mol_per_m2, self.imported_mol_per_m2, strict=True
            ):
                writer.writerow(
                    [f"{t:.6f}"] + [f"{v:.10g}" for v in totals] + [f"{v:.10g}" for v in imported]
                )
        return destination


def _digest(names: tuple[str, ...], state: NDArray[np.float64]) -> str:
    """SHA-256 over every field in a fixed byte order, independent of the machine."""
    digest = hashlib.sha256()
    for name, field in zip(names, np.asarray(state, dtype="<f8"), strict=True):
        digest.update(f"{name}:".encode())
        digest.update(np.ascontiguousarray(field).tobytes())
    return digest.hexdigest()


def time_scales(config: ReactiveTransportConfig) -> dict[str, float | None]:
    """How fast diffusion equilibrates the box against how fast anything grows.

    The ratio is recorded for information (docs/theory.md, section 5.3). The
    integrator follows transients in time, so its accuracy does not depend on
    the ratio being small, as a quasi-steady treatment's would.
    """
    domain = config.domain
    names = config.network.component_names
    height = domain.grid.size_um[-1]
    fastest = float(domain.diffusivities_um2_per_h(names).max())
    rates = [float(p.rate.maximum_per_h) for p in config.network.processes if p.rate is not None]
    if fastest <= 0 or not rates or max(rates) <= 0:
        return {"diffusion_h": None, "growth_h": None, "ratio": None}
    diffusion_h, growth_h = height**2 / fastest, 1.0 / max(rates)
    return {"diffusion_h": diffusion_h, "growth_h": growth_h, "ratio": diffusion_h / growth_h}


def build_model(config: ReactiveTransportConfig) -> ReactionTransport:
    """The discretised system a configuration describes."""
    network = config.network
    names = network.component_names
    terms = compile_rates(network)
    diffusion = Diffusion(
        config.domain.grid,
        config.domain.diffusivities_um2_per_h(names),
        config.domain.bulk(names),
    )
    return ReactionTransport(
        diffusion,
        network.stoichiometric_matrix(),
        lambda c: process_rates(terms, np.maximum(c, 0.0)),
        lambda c: rate_jacobian(terms, c),
    )


def run(
    config: ReactiveTransportConfig, *, frames: FrameWriter | None = None
) -> ReactiveTransportResult:
    """Integrate the network in space from its initial state to ``duration_h``.

    ``frames``, if given, receives (index, time_h, fields) at every recorded
    time, the initial state included. Raises
    :class:`~marse.core.simulation.ConservationError` if carbon, nitrogen or
    electrons drift.
    """
    started = datetime.now(UTC)
    network = config.network
    names = network.component_names
    grid = config.domain.grid
    model = build_model(config)
    state = config.domain.initial_state(names, config.initial_mol_per_m3, config.seed)
    ledger = Ledger(network.composition_matrix(), state, QUANTITIES)
    areal = grid.voxel_volume_um3 / grid.footprint_um2 * _AMOL_PER_UM2_TO_MOL_PER_M2
    peak = state.copy()
    imported = np.zeros(len(names))
    times, totals, imports = [0.0], [state.reshape(len(names), -1).sum(axis=1) * areal], [imported]
    if frames is not None:
        frames(0, 0.0, state)
    now, substep = 0.0, None  # the first step is estimated from the rates
    accepted = rejected = limited = failures = newton = 0
    steps = config.steps
    for step in range(1, steps + 1):
        # Times come from the step count, so they never accumulate rounding.
        target = min(step * config.timestep_h, config.duration_h)
        state, entered, stats = model.integrate(
            state,
            target - now,
            first_step=substep,
            relative_tolerance=config.relative_tolerance,
            absolute_tolerance=config.absolute_tolerance_mol_per_m3,
            peak=peak,
        )
        substep, now = stats.next_step, target
        accepted += stats.accepted
        rejected += stats.rejected
        limited += stats.limited
        failures += stats.newton_failures
        newton += stats.newton_iterations
        peak = np.maximum(peak, state)
        ledger.exchange(entered)
        ledger.check(state, step=step, time_h=now)
        imported = imported + entered
        if step % config.record_every == 0 or step == steps:
            times.append(now)
            totals.append(state.reshape(len(names), -1).sum(axis=1) * areal)
            imports.append(imported * areal)
            if frames is not None:
                frames(len(times) - 1, now, state)

    outputs: dict[str, Any] = {
        "final_time_h": now,
        "grid": {"voxels": list(grid.shape), "voxel_um": grid.voxel_um},
        "totals_mol_per_m2": {n: float(v) for n, v in zip(names, totals[-1], strict=True)},
        "imported_mol_per_m2": {n: float(v) for n, v in zip(names, imports[-1], strict=True)},
        "balance": ledger.summary(state),
        "substeps": {
            "accepted": accepted,
            "rejected": rejected,
            "limited": limited,
            "newton_failures": failures,
            "newton_iterations": newton,
        },
        "time_scales": time_scales(config),
        "final_state_sha256": _digest(names, state),
    }
    manifest = Manifest.build(
        config=config,
        models={"engine": ENGINE_VERSION, "integrator": INTEGRATOR_VERSION},
        random_streams=config.domain.random_streams,
        started=started,
        finished=datetime.now(UTC),
        steps=steps,
        outputs=outputs,
    )
    return ReactiveTransportResult(
        config=config,
        manifest=manifest,
        times_h=np.array(times),
        totals_mol_per_m2=np.array(totals),
        imported_mol_per_m2=np.array(imports),
        final_state=state,
    )
