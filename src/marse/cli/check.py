"""``marse check``: print a reaction network and whether it can run."""

from __future__ import annotations

import argparse
import math
import textwrap
from pathlib import Path

import numpy as np

from marse.chemistry import ChargeBalance
from marse.core.config import ConfigError
from marse.schemas import Network, experiment_from_dict
from marse.schemas._reading import load_json
from marse.schemas.experiment import ReactiveTransportConfig, RunConfig
from marse.schemas.formula import ELEMENT_QUANTITIES
from marse.schemas.network import RUN_FIELDS, SCHEMA_VERSION, read_document
from marse.spatial.multigrid import hierarchy


def _report_network(network: Network, source: str) -> None:
    dissolved = sum(c.phase == "dissolved" for c in network.components)
    particulate = len(network.components) - dissolved
    print(f"network     {source}, schema version {SCHEMA_VERSION}")
    if network.description:
        print(
            textwrap.fill(
                network.description, 88, initial_indent=" " * 12, subsequent_indent=" " * 12
            )
        )
    print(
        f"components  {len(network.components)}: {dissolved} dissolved, {particulate} particulate"
    )
    elements = [q for q in network.quantities if q in ELEMENT_QUANTITIES]
    symbols = "".join(f" {ELEMENT_QUANTITIES[q]:>6}" for q in elements)
    print(
        f"  {'name':<24} {'phase':<12} {'formula':<16} {'C':>6} {'N':>6} {'e-':>6}{symbols} "
        f"{'g/mol':>9}"
    )
    for c in network.components:
        carbon, nitrogen, electrons = (f"{float(q):g}" for q in c.formula.composition)
        counts = "".join(f" {float(c.formula.content(q)):>6g}" for q in elements)
        print(
            f"  {c.name:<24} {c.phase:<12} {c.formula.label():<16} {carbon:>6} {nitrogen:>6} "
            f"{electrons:>6}{counts} {c.formula.molar_mass_g_per_mol:>9.3f}"
        )
    roles = ", ".join(
        f"{c.name} {c.oxygen_role}" for c in network.components if c.oxygen_role is not None
    )
    if roles:
        print(
            textwrap.fill(
                f"oxygen roles, each checked against its processes: {roles}",
                88,
                initial_indent="  ",
                subsequent_indent="    ",
            )
        )
    if network.has_ph:
        pkw = 14.0 if network.pkw is None else float(network.pkw)
        totals = ", ".join(
            f"{c.name} ({' '.join(f'{float(k):g}' for k in c.pka)})"
            for c in network.components
            if c.pka
        )
        print(
            textwrap.fill(
                f"pH from the charge balance, pKw {pkw:g}; acid-base totals and their pKa: "
                f"{totals or 'none'}",
                88,
                initial_indent="  ",
                subsequent_indent="    ",
            )
        )
    quantities = list(network.quantities)
    conserved = ", ".join(quantities[:-1]) + " and " + quantities[-1]
    print(
        textwrap.fill(
            f"processes   {len(network.processes)}, each conserving {conserved} exactly",
            88,
            initial_indent="\n",
            subsequent_indent="            ",
        )
    )
    for p in network.processes:
        if p.growth is None:
            what = "reaction, per unit of reaction"
        else:
            what = f"growth of {p.growth.biomass} on {p.growth.substrate}, per mol formed"
        print(f"\n  {p.name} ({what})")
        print(textwrap.fill(p.equation(), 88, initial_indent="    ", subsequent_indent="      "))
        zero = [name for name in p.balanced_by if p.coefficient(name) == 0]
        if zero:
            print(f"    balanced to zero: {', '.join(zero)}")
        if p.rate is not None:
            print(
                textwrap.fill(
                    f"rate: {p.rate.describe()}",
                    88,
                    initial_indent="    ",
                    subsequent_indent="      ",
                )
            )
            if p.rate.assumed_in_excess:
                print(f"    assumed in excess: {', '.join(p.rate.assumed_in_excess)}")


def _report_domain(config: ReactiveTransportConfig) -> None:
    from marse.core.reactive_transport import build_model, time_scales

    domain = config.domain
    grid = domain.grid
    names = config.network.component_names
    size = " x ".join(f"{s:g}" for s in grid.size_um)
    voxels = " x ".join(str(n) for n in grid.shape)
    megabytes = grid.voxels * len(names) * 8 / 1e6
    print(
        f"\nspace       {grid.dimensions}-D, {voxels} voxels of {grid.voxel_um:g} um "
        f"({size} um), {grid.voxels:,} voxels, {megabytes:.3g} MB per state"
    )
    if domain.film is not None and domain.mouth is not None:
        _report_mouth(config)
    elif domain.air is None:
        bulk = ", ".join(f"{n} {v:g}" for n, v in domain.bulk_mol_per_m3.items() if v)
        print(f"  bulk liquid above, mol per m3: {bulk or 'nothing'}")
    if domain.air is not None:
        _report_air(config)
    moving = ", ".join(f"{n} {v * 1e12:g}" for n, v in domain.diffusivity_m2_per_s.items())
    print(f"  diffusivities, um2 per s: {moving or 'nothing diffuses'}")
    placed = len(domain.colonies) + sum(g.count for g in domain.random_colonies)
    kinds = sorted(
        {c.component for c in domain.colonies} | {g.component for g in domain.random_colonies}
    )
    if placed:
        print(f"  colonies: {placed} ({', '.join(kinds)})")
    if domain.plaque is not None:
        _report_plaque(config)
    if domain.surface is not None:
        _report_surface(config)
    if domain.spreading is not None and domain.spreading.mechanism == "continuum":
        _report_spreading(config)
    limit = build_model(config).diffusion.explicit_step_limit_h()
    if math.isfinite(limit):
        stability = f"an explicit step would have to be at most {limit * 3600 * 1000:.3g} ms"
    else:
        stability = "with nothing diffusing, no step is too long to be stable"
    if grid.dimensions == 1:
        print(f"  {stability}; the implicit solver solves the column directly")
    else:
        levels = hierarchy(grid.shape, len(names))
        print(
            f"  {stability}; the implicit solver uses {len(levels)} grid "
            f"level{'s' if len(levels) > 1 else ''}, "
            f"down to {' x '.join(map(str, levels[-1][0]))}"
        )
        if math.prod(levels[-1][0]) > 64:
            print(
                "  note: the coarsest level is large, which slows every solve; voxel counts "
                "divisible by 2 several times (such as 32, 48 or 64) are faster"
            )
    scales = time_scales(config)
    diffusion_h = scales["diffusion_h"]
    growth_h = scales["growth_h"]
    ratio = scales["ratio"]
    if diffusion_h is not None and growth_h is not None and ratio is not None:
        print(
            f"time scales diffusion across the box {diffusion_h * 3600:.3g} s, "
            f"fastest process {growth_h:.3g} h (ratio {ratio:.2g})"
        )


def _report_spreading(config: ReactiveTransportConfig) -> None:
    from marse.core.reactive_transport import build_model
    from marse.spatial.spreading import volume_fraction

    spreading = config.domain.spreading
    assert spreading is not None
    network = config.network
    densities = network.densities()
    rooms = ", ".join(
        f"{c.name} {float(c.density_mol_per_m3):g}"
        for c in network.components
        if c.density_mol_per_m3 is not None
    )
    assert spreading.interval_h is not None
    carried = f"; carries {', '.join(spreading.carried)}" if spreading.carried else ""
    print(
        f"  spreading: {spreading.mechanism}, every {spreading.interval_h:g} h; "
        f"packing densities, mol per m3: {rooms}{carried}"
    )
    state = config.domain.initial_state(
        network.component_names, config.initial_mol_per_m3, config.seed
    )
    phi = volume_fraction(state, densities)
    model = build_model(config)
    inverse = np.divide(1.0, densities, out=np.zeros_like(densities), where=densities > 0)
    volume_per_h = np.einsum(
        "pj,j,p...->...",
        model.stoichiometry[: len(network.processes)],
        inverse,
        model.rates(state)[: len(network.processes)],
    )
    growth = float(volume_per_h.max(initial=0.0)) * spreading.interval_h
    print(
        f"  the fullest voxel starts {float(phi.max()):.4g} full; at the starting rates, the "
        f"fastest grows by {growth:.3g} of a voxel per interval"
    )


def _report_mouth(config: ReactiveTransportConfig) -> None:
    from marse.oral import OralFluid, renewal_per_h

    domain = config.domain
    film, mouth = domain.film, domain.mouth
    assert film is not None
    assert mouth is not None
    fluid = OralFluid(mouth, film, config.network.component_names)
    rate = renewal_per_h(domain.grid, film)
    layers = round(film.thickness_um / domain.grid.voxel_um)
    mean = float(rate[..., -layers:].mean())
    print(
        f"  film        {film.thickness_um:g} um moving at {film.velocity_mm_per_min:g} mm per "
        f"minute over {film.plaque_length_mm:g} mm of plaque: renewed every "
        f"{60.0 / mean:.3g} min on average, {60.0 / rate.max():.3g} min at its surface"
    )
    cycle = (fluid.full_m3 - fluid.resting_m3) / fluid.unstimulated_m3_per_s / 60.0
    stimulus = (
        f"; up to {mouth.stimulated_flow_ml_per_min:g} more with {mouth.stimulus}"
        if (mouth.stimulus)
        else ""
    )
    if mouth.chewing_flow_ml_per_min:
        stimulus += f"; {mouth.chewing_flow_ml_per_min:g} more while chewing"
    print(
        f"  mouth       {mouth.resting_volume_ml:g} to {mouth.swallow_volume_ml:g} mL at "
        f"{mouth.unstimulated_flow_ml_per_min:g} mL per minute{stimulus}; at rest a swallow "
        f"every {cycle:.3g} min"
    )
    print(
        f"  plaque      {mouth.plaque_area_cm2:g} cm2 under {fluid.film_m3 * 1e6:.3g} mL of film; "
        f"the pool holds {fluid.thickness_um(fluid.resting_m3):.4g} um over it at rest"
    )
    saliva = ", ".join(f"{n} {v:g}" for n, v in mouth.saliva_mol_per_m3.items() if v)
    print(f"  saliva, mol per m3: {saliva or 'nothing'}")
    for intake in domain.diet:
        if intake.kind == "food":
            released = intake.released_mmol or {}
            what = ", ".join(f"{n} {v:g}" for n, v in released.items() if v)
            brings = f"releases {what} mmol" if what else "releases nothing"
            if intake.chewing:
                brings += ", chewed"
        else:
            held = intake.composition_mol_per_m3 or {}
            what = ", ".join(f"{n} {v:g}" for n, v in held.items() if v)
            brings = f"{intake.volume_ml:g} mL holding {what or 'nothing'} mol per m3"
        print(
            f"  {intake.kind:<11} at {_clock(intake.start_h)} for {intake.duration_min:g} min: "
            f"{brings}; mixes the film at {intake.mixing_per_s:g} per s"
        )
        if intake.retained is not None:
            print(
                f"              leaves {intake.retained.amount_mol_per_m2:g} mol per m2 of "
                f"{intake.retained.component} on the teeth"
            )


def _report_air(config: ReactiveTransportConfig) -> None:
    from marse.core.reactive_transport import air_of

    air = air_of(config)
    assert air is not None
    names = config.network.component_names
    held = ", ".join(
        f"{names[j]} held at {s:g} mol per m3"
        for j, s in zip(air.gases, air.saturation_mol_per_m3, strict=True)
    )
    print(f"  air         the top face is at the air: {held}")
    seconds = ", ".join(
        f"{names[j]} {3600.0 / k:.3g} s" for j, k in zip(air.gases, air.rate_per_h, strict=True)
    )
    print(f"              the top voxel follows it within: {seconds}")
    if config.domain.mouth is not None:
        print("              it holds the mouth's saliva at saturation too")


def _report_plaque(config: ReactiveTransportConfig) -> None:
    from marse.biofilm.spreading import SPREADING_VERSION

    plaque = config.domain.plaque
    assert plaque is not None
    fills = ", ".join(
        f"{c.name} at {float(c.density_mol_per_m3):g} mol per m3"
        for c in config.network.components
        if c.density_mol_per_m3 is not None
    )
    carried = f"; carries {', '.join(plaque.carried)}" if plaque.carried else ""
    print(f"  plaque      spreads up the column ({SPREADING_VERSION}): filled by {fills}{carried}")
    wear = (
        f"wears {plaque.wear_um_per_h:g} um per hour at its surface"
        if plaque.wear_um_per_h
        else "does not wear"
    )
    print(f"              at most {plaque.maximum_um:g} um high, detached above; {wear}")
    if config.domain.film is not None:
        print("              the film rides on its surface, renewed from there up")
    for cleaning in config.domain.hygiene:
        print(
            f"  {cleaning.kind:<11} at {_clock(cleaning.start_h)}: takes "
            f"{cleaning.removes_fraction:.0%} of the plaque off from its surface down"
        )


def _clock(hours: float) -> str:
    """A time in the run, in hours and minutes to a tenth of a minute."""
    whole, minutes = divmod(round(hours * 60.0, 1), 60.0)
    if not whole:
        return f"{minutes:g} min"
    return f"{whole:g} h {minutes:g} min" if minutes else f"{whole:g} h"


def _report_surface(config: ReactiveTransportConfig) -> None:
    from marse.core.reactive_transport import transfer_velocities_um_per_s
    from marse.microbes.adhesion import JAMMING_COVERAGE

    domain = config.domain
    substratum = domain.substratum
    assert substratum is not None
    assert domain.surface is not None
    assert domain.liquid is not None
    assert domain.flow is not None
    total = sum(substratum.area_um2(m) for m in substratum.materials)
    shares = ", ".join(f"{m} {substratum.area_um2(m) / total:.0%}" for m in substratum.materials)
    print(f"  substratum: {shares}; conditioning film: {domain.surface.conditioning_film}")
    print(
        f"  liquid at {domain.liquid.temperature_c:g} C, {domain.liquid.viscosity_mpa_s:g} mPa s, "
        f"sheared at {domain.flow.wall_shear_rate_per_s:g} per s, "
        f"{domain.flow.distance_from_inlet_mm:g} mm downstream"
    )
    velocities = transfer_velocities_um_per_s(config)
    for s in domain.suspension:
        arrival = velocities[s.attached] * s.cells_per_ml * 1e-12 * 1e8  # per cm2 per s
        jammed = JAMMING_COVERAGE / s.blocked_area_um2 * 1e8
        print(
            f"  {s.attached}: {s.cells_per_ml:.3g} per ml, delivered at "
            f"{velocities[s.attached]:.3g} um per s, {arrival:.3g} cells per cm2 per s; "
            f"the surface jams at {jammed:.3g} per cm2"
        )
        for m in substratum.materials:
            rule = domain.adhesion_of(s.attached, m)
            bound = rule.efficiency * arrival
            half = math.log(2) * jammed / bound / 3600 if bound > 0 else math.inf
            print(
                f"    on {m}: binds {bound:.3g} per cm2 per s at first, half-way to jamming "
                f"after {half:.3g} h without detachment or growth"
            )


def _report_ph(config: RunConfig) -> None:
    """The pH each composition a run starts from implies, from its charge balance."""
    network = config.network
    names = network.component_names
    balance = ChargeBalance.of(network)
    if isinstance(config, ReactiveTransportConfig):
        state = config.domain.initial_state(names, config.initial_mol_per_m3, config.seed)
        ph = balance.ph(state)
        mouth = config.domain.mouth
        if mouth is None:
            liquid = f"{float(balance.ph(config.domain.bulk(names))):.2f} in the bulk liquid"
        else:
            saliva = np.array([mouth.saliva_mol_per_m3.get(n, 0.0) for n in names])
            pool = np.array([mouth.initial_mol_per_m3.get(n, 0.0) for n in names])
            liquid = (
                f"{float(balance.ph(saliva)):.2f} in the saliva secreted; "
                f"{float(balance.ph(pool)):.2f} in the mouth at the start"
            )
        print(f"pH          {liquid}; {ph.min():.2f} to {ph.max():.2f} in the box at the start")
    else:
        initial = np.array([config.initial_mol_per_m3[n] for n in names])
        print(f"\npH          {float(balance.ph(initial)):.2f} at the start")


def _command_check(args: argparse.Namespace) -> int:
    path = Path(args.network)
    try:
        raw = load_json(path)
        network, values = read_document(raw)
    except ConfigError as error:
        print(f"marse: {path.name} is not a valid network: {error}")
        return 2
    _report_network(network, path.name)
    if any(key in values for key in RUN_FIELDS) or any(p.rate for p in network.processes):
        try:
            config = experiment_from_dict(raw)
        except ConfigError as error:
            print(f"\nnot yet runnable: {error}")
        else:
            if isinstance(config, ReactiveTransportConfig):
                _report_domain(config)
            if config.network.has_ph:
                _report_ph(config)
            print(
                f"\nrunnable    {config.duration_h:g} h in steps of at most "
                f"{config.timestep_h:g} h: marse run {args.network}"
            )
    return 0
