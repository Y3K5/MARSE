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

Where the domain says what its substratum is made of, cells in the liquid bind
to it as the run goes (:mod:`marse.microbes.adhesion`). The run then also
records bound cells per cm² of each material, and the area each has covered.
"""

from __future__ import annotations

import csv
import hashlib
import math
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
from marse.microbes.adhesion import AttachingSpecies, SurfaceExchange
from marse.microbes.kinetics import compile_rates, process_rates, rate_jacobian
from marse.schemas.experiment import ReactiveTransportConfig
from marse.schemas.formula import QUANTITIES
from marse.spatial.colloids import leveque_transfer_um_per_s, stokes_einstein_um2_per_s
from marse.spatial.transport import Diffusion

__all__ = ["ENGINE_VERSION", "INTEGRATOR_VERSION", "ReactiveTransportResult", "run"]

ENGINE_VERSION = "reactive_transport_v2"
INTEGRATOR_VERSION = "sdirk2_limited_v1"

type FrameWriter = Callable[[int, float, NDArray[np.float64]], None]

_AMOL_PER_UM2_TO_MOL_PER_M2 = 1e-6  # 1e-18 mol per 1e-12 m2
_PER_UM2_TO_PER_CM2 = 1e8


@dataclass(frozen=True, slots=True)
class ReactiveTransportResult:
    """A finished spatial run: areal totals over time, the final fields, the manifest."""

    config: ReactiveTransportConfig
    manifest: Manifest
    times_h: NDArray[np.float64]
    totals_mol_per_m2: NDArray[np.float64]  # (records, components)
    imported_mol_per_m2: NDArray[np.float64]  # (records, components), since the start
    final_state: NDArray[np.float64]  # (components, *voxels)
    surface: dict[str, NDArray[np.float64]] | None = None  # column -> value per record

    @property
    def component_names(self) -> tuple[str, ...]:
        return self.config.network.component_names

    def write_surface(self, path: str | Path) -> Path:
        """Write bound cells per cm² and coverage, per material, at every recorded time."""
        if self.surface is None:
            raise ValueError("this run has no substratum that cells bind to")
        destination = Path(path)
        columns = list(self.surface)
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["time_h", *columns])
            for i, t in enumerate(self.times_h):
                writer.writerow([f"{t:.6f}"] + [f"{self.surface[c][i]:.10g}" for c in columns])
        return destination

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


def transfer_velocities_um_per_s(config: ReactiveTransportConfig) -> dict[str, float]:
    """Each binding species' delivery velocity to the substratum (the Leveque flux per cell)."""
    domain = config.domain
    if domain.liquid is None or domain.flow is None:
        return {}
    velocities = {}
    for s in domain.suspension:
        diffusivity = stokes_einstein_um2_per_s(
            s.cell_diameter_um, domain.liquid.temperature_c, domain.liquid.viscosity_mpa_s
        )
        velocities[s.attached] = leveque_transfer_um_per_s(
            diffusivity, domain.flow.wall_shear_rate_per_s, domain.flow.distance_from_inlet_mm * 1e3
        )
    return velocities


def surface_exchange(config: ReactiveTransportConfig) -> SurfaceExchange | None:
    """The binding of cells to the substratum a configuration describes, if any."""
    domain = config.domain
    substratum = domain.substratum
    if substratum is None:
        return None
    names = config.network.component_names
    index = {n: i for i, n in enumerate(names)}
    velocities = transfer_velocities_um_per_s(config)
    species = []
    for s in domain.suspension:
        rules = {m: domain.adhesion_of(s.attached, m) for m in substratum.materials}
        species.append(
            AttachingSpecies(
                reversible=index[s.reversible],
                attached=index[s.attached],
                arrival_um_per_h=velocities[s.attached] * 3600.0,
                cells_per_um3=s.cells_per_ml * 1e-12,
                amol_per_cell=s.carbon_fmol_per_cell * 1e3,
                blocked_area_um2=s.blocked_area_um2,
                efficiency=substratum.per_face({m: r.efficiency for m, r in rules.items()}),
                detachment_per_h=substratum.per_face(
                    {m: r.detachment_per_h for m, r in rules.items()}
                ),
                locking_per_h=substratum.per_face({m: r.locking_per_h for m, r in rules.items()}),
            )
        )
    return SurfaceExchange(species, len(names), domain.grid.voxel_um)


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
        surface_exchange(config),
    )


def _surface_record(
    config: ReactiveTransportConfig, exchange: SurfaceExchange, state: NDArray[np.float64]
) -> dict[str, float]:
    """Bound cells per cm² of each material, species by species, and the fraction covered.

    The covered fraction is what microscopy of a surface measures: the area
    under cells, each with the footprint of its diameter, placed at random,
    1 - exp(-sum n pi d^2 / 4). It stays below one as cells pile up.
    """
    substratum = config.domain.substratum
    assert substratum is not None
    record = {}
    shade = 0.0
    for s, suspended in enumerate(config.domain.suspension):
        density = exchange.cells_per_um2(state, s)
        footprint = math.pi * suspended.cell_diameter_um**2 / 4.0
        shade = shade + density * footprint
        for material in substratum.materials:
            cells = density[substratum.mask(material)].mean() * _PER_UM2_TO_PER_CM2
            record[f"{material}_{suspended.attached}_cells_per_cm2"] = float(cells)
    covered = 1.0 - np.exp(-np.asarray(shade))
    for material in substratum.materials:
        record[f"{material}_covered"] = float(covered[substratum.mask(material)].mean())
    return record


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
    exchange = model.surface
    surface_rows = [] if exchange is None else [_surface_record(config, exchange, state)]
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
            if exchange is not None:
                surface_rows.append(_surface_record(config, exchange, state))
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
    surface = None
    if exchange is not None:
        surface = {key: np.array([row[key] for row in surface_rows]) for key in surface_rows[0]}
        outputs["surface"] = {
            "conditioning_film": config.domain.surface.conditioning_film,
            "transfer_um_per_s": transfer_velocities_um_per_s(config),
            "final": dict(surface_rows[-1]),
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
        surface=surface,
    )
