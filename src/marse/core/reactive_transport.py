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

Where the network's charges set a pH (:mod:`marse.chemistry.acid_base`), the
run records the pH over the substratum and its range in the box.

Where the domain has a salivary film and a mouth (:mod:`marse.oral`), the
film is renewed from the mouth's saliva, whose composition is solved with the
box (:mod:`marse.core.reservoir`). The run then advances in spans of at most a
minute that end at every swallow, and at the start and end of every intake of
the diet. A second ledger checks the box and the mouth together, against what
the glands secreted, what was eaten and drunk, and what was swallowed or
expelled, and the run records the mouth's volume, flow and composition.
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

from marse.chemistry import ChargeBalance
from marse.core.implicit import ReactionTransport
from marse.core.ledger import Ledger
from marse.core.provenance import Manifest
from marse.core.reservoir import ReservoirPath, ReservoirTransport
from marse.microbes.adhesion import AttachingSpecies, SurfaceExchange
from marse.microbes.kinetics import compile_rates, process_rates, rate_jacobian
from marse.oral import Diet, Inflow, OralFluid, film_layers, renewal_per_h
from marse.schemas.experiment import ReactiveTransportConfig
from marse.spatial.colloids import leveque_transfer_um_per_s, stokes_einstein_um2_per_s
from marse.spatial.transport import Diffusion

__all__ = ["ENGINE_VERSION", "INTEGRATOR_VERSION", "ReactiveTransportResult", "run"]

ENGINE_VERSION = "reactive_transport_v2"
INTEGRATOR_VERSION = "sdirk2_limited_v1"
MOUTH_VERSION = "dawes_1983_mouth_renewed_film_v1"
"""The mouth and its film, recorded in the manifest of every run that has them."""
DIET_VERSION = "rinse_drink_food_mixed_film_v1"
"""The intakes and the food they leave, recorded in the manifest of every run with a diet."""

type FrameWriter = Callable[[int, float, NDArray[np.float64]], None]

_AMOL_PER_UM2_TO_MOL_PER_M2 = 1e-6  # 1e-18 mol per 1e-12 m2
_PER_UM2_TO_PER_CM2 = 1e8
_LONGEST_SPAN_S = 60.0  # the mouth runs ahead of the box at most this far
_ML_PER_M3 = 1e6


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
    ph: dict[str, NDArray[np.float64]] | None = None  # column -> value per record
    mouth: dict[str, NDArray[np.float64]] | None = None  # column -> value per record

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

    def write_ph(self, path: str | Path) -> Path:
        """Write the pH at the substratum and its range in the box, at every recorded time."""
        if self.ph is None:
            raise ValueError("this run's network sets no pH")
        destination = Path(path)
        columns = list(self.ph)
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["time_h", *columns])
            for i, t in enumerate(self.times_h):
                writer.writerow([f"{t:.6f}"] + [f"{self.ph[c][i]:.6f}" for c in columns])
        return destination

    def write_mouth(self, path: str | Path) -> Path:
        """Write the mouth's volume, flow, swallows and composition at every recorded time."""
        if self.mouth is None:
            raise ValueError("this run has no mouth")
        destination = Path(path)
        columns = list(self.mouth)
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["time_h", *columns])
            for i, t in enumerate(self.times_h):
                writer.writerow([f"{t:.6f}"] + [f"{self.mouth[c][i]:.10g}" for c in columns])
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


def _digest(
    names: tuple[str, ...], state: NDArray[np.float64], pool: NDArray[np.float64] | None = None
) -> str:
    """SHA-256 over every field in a fixed byte order, independent of the machine.

    A run with a mouth adds the mouth's composition, so a replay reproduces both.
    """
    digest = hashlib.sha256()
    for name, field in zip(names, np.asarray(state, dtype="<f8"), strict=True):
        digest.update(f"{name}:".encode())
        digest.update(np.ascontiguousarray(field).tobytes())
    if pool is not None:
        digest.update(b"mouth:")
        digest.update(np.ascontiguousarray(np.asarray(pool, dtype="<f8")).tobytes())
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
    """The discretised system a configuration describes.

    Under a film, a :class:`~marse.core.reservoir.ReservoirTransport`, closed at
    the top and bordered by the mouth's pool.
    """
    network = config.network
    names = network.component_names
    terms = compile_rates(network)
    domain = config.domain

    def rates(c: NDArray[np.float64]) -> NDArray[np.float64]:
        return process_rates(terms, np.maximum(c, 0.0))

    def jacobian(c: NDArray[np.float64]) -> NDArray[np.float64]:
        return rate_jacobian(terms, c)

    if domain.film is not None and domain.mouth is not None:
        fluid = OralFluid(domain.mouth, domain.film, names)
        diffusion = Diffusion(
            domain.grid,
            domain.diffusivities_um2_per_h(names),
            np.zeros(len(names)),
            closed_top=True,
        )
        return ReservoirTransport(
            diffusion,
            network.stoichiometric_matrix(),
            rates,
            jacobian,
            exchanged=[j for j, c in enumerate(network.components) if c.phase == "dissolved"],
            reference_um=fluid.reference_um,
            surface=surface_exchange(config),
        )
    diffusion = Diffusion(
        domain.grid,
        domain.diffusivities_um2_per_h(names),
        domain.bulk(names),
    )
    return ReactionTransport(
        diffusion,
        network.stoichiometric_matrix(),
        rates,
        jacobian,
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


def _ph_record(balance: ChargeBalance, state: NDArray[np.float64]) -> dict[str, float]:
    """The pH over the substratum, the bottom layer of voxels, and its range in the whole box."""
    ph = balance.ph(state)
    bottom = ph[..., 0]
    return {
        "substratum_mean": float(bottom.mean()),
        "substratum_min": float(bottom.min()),
        "substratum_max": float(bottom.max()),
        "box_min": float(ph.min()),
        "box_max": float(ph.max()),
    }


def run(
    config: ReactiveTransportConfig, *, frames: FrameWriter | None = None
) -> ReactiveTransportResult:
    """Integrate the network in space from its initial state to ``duration_h``.

    ``frames``, if given, receives (index, time_h, fields) at every recorded
    time, the initial state included. Raises
    :class:`~marse.core.simulation.ConservationError` if carbon, nitrogen or
    electrons drift.
    """
    if config.domain.mouth is not None:
        return _run_with_mouth(config, frames)
    started = datetime.now(UTC)
    network = config.network
    names = network.component_names
    grid = config.domain.grid
    model = build_model(config)
    state = config.domain.initial_state(names, config.initial_mol_per_m3, config.seed)
    ledger = Ledger(network.composition_matrix(), state, network.quantities)
    balance = ChargeBalance.of(network) if network.has_ph else None
    ph_rows = [] if balance is None else [_ph_record(balance, state)]
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
            absolute_tolerance=config.absolute_tolerances(),
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
            if balance is not None:
                ph_rows.append(_ph_record(balance, state))
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
    ph = None
    if balance is not None:
        ph = {key: np.array([row[key] for row in ph_rows]) for key in ph_rows[0]}
        lowest = int(np.argmin(ph["substratum_min"]))
        outputs["ph"] = {
            "final": dict(ph_rows[-1]),
            "lowest_at_substratum": float(ph["substratum_min"][lowest]),
            "lowest_at_substratum_h": float(times[lowest]),
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
        ph=ph,
    )


def _path(
    fluid: OralFluid,
    start_h: float,
    stretch: Any,
    secreted: NDArray[np.float64],
    inflow: Inflow,
) -> ReservoirPath:
    """The pool's thickness over one span, from the mouth's run ahead, and what enters it."""
    supply = None
    if inflow.released_mol_per_s is not None:  # mol/s into the mouth; mol/m3 x um/h per area
        supply = inflow.released_mol_per_s * 3600.0 / fluid.area_m2 * 1e6
    return ReservoirPath(
        start_h=start_h,
        times_h=stretch.times_s / 3600.0,
        thickness_um=np.array([fluid.thickness_um(v) for v in stretch.volumes_m3]),
        growth_um_per_h=np.array([fluid.growth_um_per_h(q) for q in stretch.flows_m3_per_s]),
        secreted_mol_per_m3=secreted,
        drink_um_per_h=fluid.growth_um_per_h(inflow.liquid_m3_per_s),
        drink_mol_per_m3=inflow.liquid_mol_per_m3,
        supply_per_h=supply,
    )


def _mouth_record(
    fluid: OralFluid,
    engine: ReservoirTransport,
    pool: NDArray[np.float64],
    names: tuple[str, ...],
    balance: ChargeBalance | None,
) -> dict[str, float]:
    """The mouth's volume, flow and swallows, and the composition of its pool."""
    concentration = pool * engine.reference_um / fluid.thickness_um(fluid.volume_m3)
    stimulus = 0.0 if fluid.stimulus is None else float(concentration[fluid.stimulus])
    record = {
        "volume_ml": fluid.volume_m3 * _ML_PER_M3,
        "flow_ml_per_min": fluid.flow_m3_per_s(stimulus) * _ML_PER_M3 * 60.0,
        "swallows": float(fluid.swallows),
    }
    for j in engine.exchanged:
        record[f"{names[j]}_mol_per_m3"] = float(concentration[j])
    if balance is not None:
        record["ph"] = float(balance.ph(concentration))
    return record


def _run_with_mouth(
    config: ReactiveTransportConfig, frames: FrameWriter | None
) -> ReactiveTransportResult:
    """Integrate a box under a salivary film, bordered by the mouth's saliva.

    Each span runs the mouth ahead to its end, at most a minute, or to the
    next swallow; the box and the pool are then integrated together over it,
    and a swallow takes the pool back to its resting volume. Spans also end
    where an intake starts or ends: a rinse is taken in, or expelled, and food
    left on the teeth is placed in the film. The box's ledger books what
    crossed into the film and the food placed in it; the second ledger checks
    the box and the mouth together, against what was secreted, eaten,
    swallowed and expelled.
    """
    started = datetime.now(UTC)
    network = config.network
    names = network.component_names
    domain = config.domain
    assert domain.film is not None
    assert domain.mouth is not None
    grid = domain.grid
    engine = build_model(config)
    assert isinstance(engine, ReservoirTransport)
    fluid = OralFluid(domain.mouth, domain.film, names)
    renewal = renewal_per_h(grid, domain.film)
    in_film = np.zeros(grid.shape)
    in_film[..., -film_layers(grid, domain.film) :] = 1.0
    diet = Diet(domain.diet, names)
    box = domain.initial_state(names, config.initial_mol_per_m3, config.seed)
    pool = np.array([domain.mouth.initial_mol_per_m3.get(n, 0.0) for n in names])
    y = engine.pack(box, pool)  # at the resting volume, the pool's unknowns are its concentrations
    areal = grid.voxel_volume_um3 / grid.footprint_um2 * _AMOL_PER_UM2_TO_MOL_PER_M2
    to_mol_per_m2 = engine.reference_um * 1e-6  # the pool's unknowns, per unit area

    def everything(box: NDArray[np.float64], pool: NDArray[np.float64]) -> NDArray:
        return box.reshape(len(names), -1).sum(axis=1) * areal + pool * to_mol_per_m2

    ledger = Ledger(network.composition_matrix(), box, network.quantities)
    whole = Ledger(network.composition_matrix(), everything(box, pool), network.quantities)
    balance = ChargeBalance.of(network) if network.has_ph else None
    peak = y.copy()
    imported = np.zeros(len(names))
    secreted = np.zeros(len(names))
    swallowed = np.zeros(len(names))
    eaten = np.zeros(len(names))
    expelled = np.zeros(len(names))
    times, totals, imports = [0.0], [box.reshape(len(names), -1).sum(axis=1) * areal], [imported]
    exchange = engine.surface
    surface_rows = [] if exchange is None else [_surface_record(config, exchange, box)]
    ph_rows = [] if balance is None else [_ph_record(balance, box)]
    mouth_rows = [_mouth_record(fluid, engine, pool, names, balance)]
    if frames is not None:
        frames(0, 0.0, box)
    now, substep = 0.0, None
    accepted = rejected = limited = failures = newton = spans = 0
    returned = 0.0  # how fast the box gave the stimulus back to the pool over the last span, mol/s
    steps = config.steps
    for step in range(1, steps + 1):
        target = min(step * config.timestep_h, config.duration_h)
        while target - now > 1e-12 * target:
            # An intake starts or ends: a rinse is taken in or expelled, and food left on
            # the teeth is placed in the film.
            events = diet.due(now)
            if events:
                returned = 0.0  # what the box did before an intake says nothing about after it
                box, pool = engine.unpack(y)
                for event in events:
                    intake = event.intake
                    if event.starts and intake.kind == "rinse":
                        volume, amounts = diet.rinse(intake)
                        fluid.take(volume)
                        brought = amounts / fluid.area_m2  # mol/m2
                        eaten += brought
                        whole.exchange(brought)
                        pool = pool + brought / to_mol_per_m2
                    if not event.starts and intake.kind == "rinse":
                        kept = fluid.expel()
                        gone = pool * (1.0 - kept) * to_mol_per_m2
                        expelled += gone
                        whole.exchange(-gone)
                        pool = pool * kept
                    pocket = None if event.starts else diet.pocket(intake, grid, domain.film)
                    if pocket is not None:
                        box = box + pocket
                        placed = pocket.reshape(len(names), -1).sum(axis=1)
                        ledger.exchange(placed)
                        imported = imported + placed
                        eaten += placed * areal
                        whole.exchange(placed * areal)
                y = engine.pack(box, pool)
                ledger.check(box, step=step, time_h=now)
                whole.check(everything(box, pool), step=step, time_h=now)
                peak = np.maximum(peak, y)
            inflow = diet.inflow()
            engine.exchange_per_h = (
                renewal if inflow.mixing_per_h == 0.0 else renewal + inflow.mixing_per_h * in_film
            )
            box, pool = engine.unpack(y)
            stimulus = 0.0
            if fluid.stimulus is not None:
                stimulus = pool[fluid.stimulus] * to_mol_per_m2 * fluid.area_m2
            stop = min(target, diet.next_h())
            left = (stop - now) * 3600.0
            stretch = fluid.run_ahead(stimulus, min(_LONGEST_SPAN_S, left), inflow, returned)
            # Land exactly on the step's end, or the intake's, when the span reaches it, so
            # times never drift.
            end = stop if stretch.seconds >= left * (1 - 1e-12) else now + stretch.seconds / 3600.0
            span = end - now
            if span > 0:
                growth = stretch.volumes_m3[-1] - stretch.volumes_m3[0]
                saliva = growth / stretch.seconds - inflow.liquid_m3_per_s
                # The span keeps its own clock, from zero. On the run's clock, a step of
                # seconds an hour or a day in would be the difference of two large times,
                # which loses digits in proportion to the time: enough, over a day of
                # meals, to leave 4e-13 of the sugar eaten unaccounted for.
                engine.path = _path(fluid, 0.0, stretch, fluid.secreted(saliva), inflow)
                y, entered, stats = engine.integrate(
                    y,
                    span,
                    first_step=substep,
                    relative_tolerance=config.relative_tolerance,
                    absolute_tolerance=config.absolute_tolerances(),
                    peak=peak,
                    start_h=0.0,
                )
                substep = stats.next_step
                accepted += stats.accepted
                rejected += stats.rejected
                limited += stats.limited
                failures += stats.newton_failures
                newton += stats.newton_iterations
                spans += 1
                ledger.exchange(entered)
                imported = imported + entered
                if fluid.stimulus is not None:
                    given = entered[fluid.stimulus] * areal * fluid.area_m2  # mol, into the box
                    returned = -given / (span * 3600.0)
                added = engine.path.added(0.0, span) * 1e-6  # mol/m2
                taken = engine.path.taken(0.0, span) * 1e-6
                secreted += added - taken
                eaten += taken
                whole.exchange(added)
            kept = fluid.end(stretch)
            box, pool = engine.unpack(y)
            if kept < 1.0:
                gone = pool * (1.0 - kept) * to_mol_per_m2
                swallowed += gone
                whole.exchange(-gone)
                pool = pool * kept
                y = engine.pack(box, pool)
            ledger.check(box, step=step, time_h=end)
            whole.check(everything(box, pool), step=step, time_h=end)
            peak = np.maximum(peak, y)
            now = end
        now = target
        if step % config.record_every == 0 or step == steps:
            box, pool = engine.unpack(y)
            times.append(now)
            totals.append(box.reshape(len(names), -1).sum(axis=1) * areal)
            imports.append(imported * areal)
            if exchange is not None:
                surface_rows.append(_surface_record(config, exchange, box))
            if balance is not None:
                ph_rows.append(_ph_record(balance, box))
            mouth_rows.append(_mouth_record(fluid, engine, pool, names, balance))
            if frames is not None:
                frames(len(times) - 1, now, box)

    box, pool = engine.unpack(y)
    concentration = pool * engine.reference_um / fluid.thickness_um(fluid.volume_m3)

    def per_component(amounts: NDArray[np.float64]) -> dict[str, float]:
        return {n: float(v) for n, v in zip(names, amounts, strict=True)}

    mouth: dict[str, Any] = {
        "spans": spans,
        "swallows": fluid.swallows,
        "final_volume_ml": fluid.volume_m3 * _ML_PER_M3,
        "final_mol_per_m3": {names[j]: float(concentration[j]) for j in engine.exchanged},
        "secreted_mol_per_m2": per_component(secreted),
        "swallowed_mol_per_m2": per_component(swallowed),
    }
    models = {"engine": ENGINE_VERSION, "integrator": INTEGRATOR_VERSION, "mouth": MOUTH_VERSION}
    if domain.diet:
        models["diet"] = DIET_VERSION
        mouth["intakes"] = diet.taken
        mouth["eaten_mol_per_m2"] = per_component(eaten)
        mouth["expelled_mol_per_m2"] = per_component(expelled)
    mouth["balance"] = whole.summary(everything(box, pool))
    outputs: dict[str, Any] = {
        "final_time_h": now,
        "grid": {"voxels": list(grid.shape), "voxel_um": grid.voxel_um},
        "totals_mol_per_m2": {n: float(v) for n, v in zip(names, totals[-1], strict=True)},
        "imported_mol_per_m2": {n: float(v) for n, v in zip(names, imports[-1], strict=True)},
        "balance": ledger.summary(box),
        "substeps": {
            "accepted": accepted,
            "rejected": rejected,
            "limited": limited,
            "newton_failures": failures,
            "newton_iterations": newton,
        },
        "time_scales": time_scales(config),
        "mouth": mouth,
        "final_state_sha256": _digest(names, box, pool),
    }
    surface = None
    if exchange is not None:
        surface = {key: np.array([row[key] for row in surface_rows]) for key in surface_rows[0]}
        outputs["surface"] = {
            "conditioning_film": domain.surface.conditioning_film if domain.surface else "",
            "transfer_um_per_s": transfer_velocities_um_per_s(config),
            "final": dict(surface_rows[-1]),
        }
    ph = None
    if balance is not None:
        ph = {key: np.array([row[key] for row in ph_rows]) for key in ph_rows[0]}
        lowest = int(np.argmin(ph["substratum_min"]))
        outputs["ph"] = {
            "final": dict(ph_rows[-1]),
            "lowest_at_substratum": float(ph["substratum_min"][lowest]),
            "lowest_at_substratum_h": float(times[lowest]),
        }
    manifest = Manifest.build(
        config=config,
        models=models,
        random_streams=domain.random_streams,
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
        final_state=box,
        surface=surface,
        ph=ph,
        mouth={key: np.array([row[key] for row in mouth_rows]) for key in mouth_rows[0]},
    )
