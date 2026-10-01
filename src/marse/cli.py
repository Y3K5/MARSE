"""Command-line entry point: ``marse``.

Four commands matter at this stage:

    marse run experiment.json        run a batch or biofilm simulation, or a
                                     version 2 reaction network, in a closed
                                     box or, with a domain, in space
    marse ecosystem experiment.json  run a 2-D multispecies ecosystem
    marse replay manifest.json       re-run any of them from its manifest and compare
    marse check network.json         check a reaction network (schema version 2)

``replay`` is the reproducibility claim made executable: it rebuilds the
configuration from the manifest, verifies the checksum, runs again, and reports
whether the results match.
"""

from __future__ import annotations

import argparse
import json
import math
import textwrap
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from marse import __version__
from marse.analysis.ensemble import ScenarioBatch, run_batch
from marse.chemistry import ChargeBalance
from marse.core.config import ConfigError, load_experiment
from marse.core.framestore import FrameStore as FieldStore
from marse.core.provenance import Manifest
from marse.core.reactive_transport import ReactiveTransportResult
from marse.core.reactive_transport import run as run_reactive_transport
from marse.core.simulation import BiofilmProfileResult, SimulationResult, run
from marse.core.well_mixed import WellMixedResult
from marse.core.well_mixed import run as run_well_mixed
from marse.ecosystem import EcosystemConfig, EcosystemResult, write_viewer
from marse.ecosystem import load_experiment as load_ecosystem_experiment
from marse.ecosystem import run as run_ecosystem
from marse.ecosystem.framestore import EcosystemFrameSink, FrameStore, StoredFrames
from marse.ecosystem.model import recorded_steps
from marse.microbes.niche import NicheError, load_niche_scan, run_niche_scan
from marse.schemas import Network, experiment_from_dict
from marse.schemas._reading import load_json
from marse.schemas.experiment import ReactiveTransportConfig, RunConfig
from marse.schemas.experiment import load_experiment as load_network_experiment
from marse.schemas.formula import ELEMENT_QUANTITIES
from marse.schemas.network import RUN_FIELDS, SCHEMA_VERSION, read_document
from marse.spatial.multigrid import hierarchy
from marse.spatial.vtk import write_pvd, write_vti

Result = SimulationResult | BiofilmProfileResult | WellMixedResult | ReactiveTransportResult


def _write_outputs(result: Result, output_dir: Path) -> tuple[Path, Path]:
    """Write the run's data file and its manifest. The data file differs by kind."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if isinstance(result, BiofilmProfileResult):
        data = result.write_profile(output_dir / "profile.csv")
    elif isinstance(result, ReactiveTransportResult):
        data = result.write_totals(output_dir / "totals.csv")
    else:
        data = result.write_trajectory(output_dir / "trajectory.csv")
    manifest = result.manifest.write(output_dir / "manifest.json")
    return data, manifest


def _summarise(result: Result) -> None:
    print(f"run_id      {result.manifest.run_id}")
    print(f"experiment  {result.config.experiment_id}")
    print(f"kind        {result.config.kind}")

    if isinstance(result, BiofilmProfileResult):
        outputs = result.manifest.outputs
        settings = result.config.biofilm
        print(f"biofilm     {settings.thickness_um:g} um in {settings.cells} nodes")
        print(f"solved      {outputs['newton_iterations']} Newton iterations")
        print(f"penetration {outputs['penetration_depth_um']:.1f} um")
        print(f"active zone {outputs['active_zone_um']:.1f} um")
        for name, rate in outputs["mean_growth_rate_per_h"].items():
            surface = outputs["surface_growth_rate_per_h"][name]
            share = rate / surface if surface else 0.0
            print(f"  {name:<24} mean {rate:.6g} /h ({share:.1%} of the surface rate)")
        return

    final = result.final_state
    print(f"steps       {final.step} over {final.time_h:g} h")
    print(
        f"substrate   {final.substrate_mm:.6g} mM left of {result.config.substrate.initial_mm:g} mM"
    )
    for name, value in zip(result.organism_names, final.biomass_g_per_l, strict=True):
        print(f"  {name:<24} {value:.6g} g/L")


_BIOFILM_KEYS = ("penetration_depth_um", "active_zone_um", "base_concentration_mm")


def _summarise_well_mixed(result: WellMixedResult) -> None:
    outputs = result.manifest.outputs
    substeps = outputs["substeps"]
    print(f"run_id      {result.manifest.run_id}")
    print(f"experiment  {result.config.experiment_id}")
    print(f"kind        {result.config.kind}")
    print(
        f"steps       {result.config.steps} over {outputs['final_time_h']:g} h: "
        f"{substeps['accepted']} adaptive substeps, {substeps['rejected']} retried, "
        f"{substeps['limited']} limited to keep concentrations positive"
    )
    print("final       mol per m3")
    for name, value in outputs["final_mol_per_m3"].items():
        print(f"  {name:<24} {value:.6g}")
    if "ph" in outputs:
        ph = outputs["ph"]
        print(
            f"pH          {ph['final']:.3f} at the end; lowest {ph['lowest']:.3f} "
            f"at {ph['lowest_at_h']:g} h"
        )
    worst = max(q["largest_relative_residual"] for q in outputs["balance"].values())
    print(f"balance     {_listed(outputs['balance'])} conserved to {worst:.1e} of their totals")


def _listed(quantities: Sequence[str]) -> str:
    """``carbon, nitrogen and electrons``, from the quantities a balance lists."""
    names = list(quantities)
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]


def _summarise_in_space(result: ReactiveTransportResult) -> None:
    outputs = result.manifest.outputs
    substeps = outputs["substeps"]
    grid = result.config.domain.grid
    size = " x ".join(f"{s:g}" for s in grid.size_um)
    voxels = " x ".join(str(n) for n in grid.shape)
    print(f"run_id      {result.manifest.run_id}")
    print(f"experiment  {result.config.experiment_id}")
    print(f"kind        {result.config.kind}")
    print(f"space       {grid.dimensions}-D, {voxels} voxels of {grid.voxel_um:g} um ({size} um)")
    print(
        f"steps       {result.config.steps} over {outputs['final_time_h']:g} h: "
        f"{substeps['accepted']} implicit substeps, {substeps['rejected']} retried, "
        f"{substeps['limited']} limited to keep concentrations positive"
    )
    faces = ["the film" if "mouth" in outputs else "the air" if "air" in outputs else "the top"]
    if "mouth" in outputs and "air" in outputs:
        faces.append("the air")
    if "surface" in outputs:
        faces.append("the substratum")
    boundary = _listed(faces)
    if "mouth" in outputs and any(intake.retained for intake in result.config.domain.diet):
        boundary += " (food left on the teeth included)"
    print(f"{'final':<24}   mol per m2 of surface   entered through {boundary}")
    for name, value in outputs["totals_mol_per_m2"].items():
        entered = outputs["imported_mol_per_m2"][name]
        print(f"  {name:<24} {value:>12.6g}   {entered:>+12.4g}")
    if "surface" in outputs:
        _summarise_surface(result)
    if "mouth" in outputs:
        _summarise_mouth(result)
    if "plaque" in outputs:
        _summarise_plaque(result)
    if "air" in outputs:
        _summarise_air(result)
    if "ph" in outputs:
        ph = outputs["ph"]
        final = ph["final"]
        print(
            f"pH          at the substratum {final['substratum_mean']:.3f} at the end "
            f"(lowest {ph['lowest_at_substratum']:.3f} at {ph['lowest_at_substratum_h']:g} h); "
            f"in the box {final['box_min']:.3f} to {final['box_max']:.3f}"
        )
    worst = max(q["largest_relative_residual"] for q in outputs["balance"].values())
    print(
        f"balance     {_listed(outputs['balance'])}, counting what crossed {boundary}, "
        f"conserved to {worst:.1e}"
    )


def _summarise_mouth(result: ReactiveTransportResult) -> None:
    mouth = result.manifest.outputs["mouth"]
    worst = max(q["largest_relative_residual"] for q in mouth["balance"].values())
    print(
        f"mouth       {mouth['swallows']} swallows in {mouth['spans']} spans; "
        f"{mouth['final_volume_ml']:.3g} mL at the end"
    )
    print(f"  {'saliva at the end':<24} mol per m3")
    for name, value in mouth["final_mol_per_m3"].items():
        print(f"  {name:<24} {value:>12.6g}")
    if result.mouth is not None and "ph" in result.mouth:
        print(f"  {'pH':<24} {result.mouth['ph'][-1]:>12.3f}")
    counted = ["secreted", "swallowed"]
    if "intakes" in mouth:
        counted = ["secreted", "eaten", "swallowed"]
        print(f"  {'taken in':<24} mol per m2 of plaque, from {mouth['intakes']} intake(s)")
        for name, value in mouth["eaten_mol_per_m2"].items():
            if value:
                print(f"  {name:<24} {value:>12.6g}")
    if "expelled_mol_per_m2" in mouth:
        counted.append("expelled")
    if "air" in result.manifest.outputs:
        counted.append("exchanged with the air")
    print(
        f"  box and mouth together, counting what was {_listed(counted)}, conserved to {worst:.1e}"
    )


def _summarise_plaque(result: ReactiveTransportResult) -> None:
    plaque = result.manifest.outputs["plaque"]
    assert result.plaque is not None
    thickness = result.plaque["thickness_um"]
    print(
        f"plaque      {plaque['final_thickness_um']:.4g} um at the end, "
        f"{thickness.min():.4g} to {thickness.max():.4g} um over the run "
        f"(spreading {plaque['spreading']})"
    )
    detached = ", ".join(f"{n} {v:.4g}" for n, v in plaque["detached_mol_per_m2"].items() if v)
    print(f"  detached, worn and over the top, mol per m2: {detached or 'nothing'}")
    if "removed_mol_per_m2" in plaque:
        cleanings = result.manifest.outputs["mouth"].get("cleanings", 0)
        removed = ", ".join(f"{n} {v:.4g}" for n, v in plaque["removed_mol_per_m2"].items() if v)
        print(f"  taken off by {cleanings} cleaning(s), mol per m2: {removed or 'nothing'}")


def _summarise_air(result: ReactiveTransportResult) -> None:
    air = result.manifest.outputs["air"]
    gave = ", ".join(f"{n} {v:.4g}" for n, v in air["exchanged_mol_per_m2"].items())
    into = "the box and the mouth" if "mouth" in result.manifest.outputs else "the box"
    print(f"air         gave {gave} mol per m2 to {into}, less what it took back")


def _summarise_surface(result: ReactiveTransportResult) -> None:
    domain = result.config.domain
    final = result.manifest.outputs["surface"]["final"]
    substratum = domain.substratum
    assert substratum is not None
    species = [s.attached for s in domain.suspension]
    print(f"surface     {domain.surface.conditioning_film if domain.surface else ''}")
    print(f"  {'bound cells per cm2':<22}" + "".join(f"{s:>14}" for s in species) + "   covered")
    for material in substratum.materials:
        cells = "".join(f"{final[f'{material}_{s}_cells_per_cm2']:>14.4g}" for s in species)
        print(f"  {material:<22}{cells}   {final[f'{material}_covered']:>6.1%}")


def _comparable(outputs: dict) -> dict[str, float | str]:
    """The values a replay must reproduce, flattened so every kind compares alike."""
    if "totals_mol_per_m2" in outputs:  # a reaction network in space
        values: dict[str, float | str] = {
            "final_time_h": float(outputs["final_time_h"]),
            "final state digest": str(outputs["final_state_sha256"]),
        }
        for group in ("totals_mol_per_m2", "imported_mol_per_m2"):
            values.update({f"{group}[{n}]": float(v) for n, v in outputs[group].items()})
        if "surface" in outputs:
            final = outputs["surface"]["final"]
            values.update({f"surface[{k}]": float(v) for k, v in final.items()})
        return values
    if "final_mol_per_m3" in outputs:  # a reaction network in a closed box
        values = {
            "final_time_h": float(outputs["final_time_h"]),
            "final state digest": str(outputs["final_state_sha256"]),
        }
        values.update({f"final[{n}]": float(v) for n, v in outputs["final_mol_per_m3"].items()})
        return values
    if "final_state_sha256" in outputs:  # an ecosystem run
        values = {
            "final_time_h": float(outputs["final_time_h"]),
            "final state digest": str(outputs["final_state_sha256"]),
        }
        for group in ("biomass_summed_over_cells", "nutrient_summed_over_cells"):
            values.update({f"{group}[{n}]": float(v) for n, v in outputs[group].items()})
        values.update(
            {f"mutated_cells[{n}]": float(v) for n, v in outputs["mutated_cells"].items()}
        )
        return values
    if "final_state" in outputs:  # a batch run
        final = outputs["final_state"]
        values = {"substrate": float(final["substrate_mm"])}
        values.update({name: float(v) for name, v in final["biomass_g_per_l"].items()})
        return values
    values = {key: float(outputs[key]) for key in _BIOFILM_KEYS}
    values.update(
        {f"mean growth of {n}": float(v) for n, v in outputs["mean_growth_rate_per_h"].items()}
    )
    return values


def _is_version_2(path: str) -> bool:
    """Whether a configuration file declares a schema version (only version 2 does)."""
    try:
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False  # the regular loader reports the problem
    return isinstance(raw, dict) and "schema_version" in raw


def _run_in_space(
    config: ReactiveTransportConfig, output_dir: Path
) -> tuple[ReactiveTransportResult, tuple[Path, ...]]:
    """Run in space, streaming every recorded frame to a frame store and to VTK files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    grid = config.domain.grid
    names = config.network.component_names
    records = 1 + sum(
        1
        for step in range(1, config.steps + 1)
        if step % config.record_every == 0 or step == config.steps
    )
    balance = ChargeBalance.of(config.network) if config.network.has_ph else None
    layout = {f"{n}_mol_per_m3": (grid.shape, "<f4") for n in names}
    if balance is not None:
        layout["ph"] = (grid.shape, "<f4")
    store = FieldStore.create(
        output_dir / "frames",
        capacity=records,
        fields=layout,
        metadata={"voxel_um": grid.voxel_um, "axes": ["x", "y", "z"][-grid.dimensions :]},
        overwrite=True,
    )
    vtk = output_dir / "vtk"
    vtk.mkdir(exist_ok=True)
    series: list[tuple[float, str]] = []

    def record(index: int, time_h: float, state: np.ndarray) -> None:
        fields = {f"{n}_mol_per_m3": state[j] for j, n in enumerate(names)}
        if balance is not None:
            fields["ph"] = balance.ph(state)
        store.append(time_h=time_h, step=index, arrays=fields)
        name = f"frame_{index:04d}.vti"
        write_vti(vtk / name, grid, fields)
        series.append((time_h, name))

    try:
        result = run_reactive_transport(config, frames=record)
    finally:
        index = store.close()
    pvd = write_pvd(vtk / "run.pvd", series)
    totals, manifest = _write_outputs(result, output_dir)
    written = [totals]
    if result.surface is not None:
        written.append(result.write_surface(output_dir / "surface.csv"))
    if result.ph is not None:
        written.append(result.write_ph(output_dir / "ph.csv"))
    if result.mouth is not None:
        written.append(result.write_mouth(output_dir / "mouth.csv"))
    if result.plaque is not None:
        written.append(result.write_plaque(output_dir / "plaque.csv"))
    return result, (*written, index, pvd, manifest)


def _run_network(args: argparse.Namespace) -> int:
    config = load_network_experiment(args.experiment)
    if isinstance(config, ReactiveTransportConfig):
        output_dir = (
            Path(args.output)
            if args.output
            else Path(args.experiment).parent / "runs" / Manifest.run_id_for(config)
        )
        result_in_space, written = _run_in_space(config, output_dir)
        _summarise_in_space(result_in_space)
        print()
        for path in written:
            print(f"wrote {path}")
        pvd = next(path for path in written if path.suffix == ".pvd")
        print(f"\nopen {pvd} in ParaView to explore the run in 3-D")
        print(f"replay it with:  marse replay {written[-1]}")
        return 0
    result = run_well_mixed(config)
    output_dir = (
        Path(args.output)
        if args.output
        else Path(args.experiment).parent / "runs" / result.manifest.run_id
    )
    trajectory, manifest = _write_outputs(result, output_dir)
    _summarise_well_mixed(result)
    print(f"\nwrote {trajectory}")
    print(f"wrote {manifest}")
    print(f"\nreplay it with:  marse replay {manifest}")
    return 0


def _command_run(args: argparse.Namespace) -> int:
    if _is_version_2(args.experiment):
        return _run_network(args)
    config = load_experiment(args.experiment)
    result = run(config)
    output_dir = (
        Path(args.output)
        if args.output
        else Path(args.experiment).parent / "runs" / result.manifest.run_id
    )
    trajectory, manifest = _write_outputs(result, output_dir)
    _summarise(result)
    print(f"\nwrote {trajectory}")
    print(f"wrote {manifest}")
    print(f"\nreplay it with:  marse replay {manifest}")
    return 0


MAX_STORED_FRAMES = 200
"""About how many frames ``marse ecosystem`` stores by default, whatever the run length."""


def _frame_every(config: EcosystemConfig, interval_h: float | None) -> int:
    """Steps between stored frames: from ``--frame-interval-h``, else a cap on the count."""
    if interval_h is not None:
        if not interval_h > 0:
            raise ValueError("--frame-interval-h must be positive")
        return max(1, round(interval_h / config.timestep_h))
    return max(1, math.ceil(config.steps / MAX_STORED_FRAMES))


def _run_ecosystem_into(
    config: EcosystemConfig, output_dir: Path, frame_every: int
) -> tuple[EcosystemResult, int, tuple[Path, ...]]:
    """Run, streaming frames into ``output_dir/frames``; then write the viewer and manifest.

    Memory stays flat however long the run: frames go to disk as they are made.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    capacity = len(recorded_steps(config.steps, frame_every))
    sink = EcosystemFrameSink.create(
        output_dir / "frames", config, capacity=capacity, overwrite=True
    )
    try:
        result = run_ecosystem(config, frame_every=frame_every, sink=sink)
    finally:
        index = sink.close()
    store = FrameStore.open(output_dir / "frames")
    try:
        viewer = write_viewer(result, output_dir / "viewer.html", frames=StoredFrames(store))
    finally:
        store.close()  # release the memory maps, which Windows requires before files move
    manifest = result.manifest.write(output_dir / "manifest.json")
    return result, capacity, (index, viewer, manifest)


def _command_replay(args: argparse.Namespace) -> int:
    original = Manifest.read(args.manifest)
    config = original.experiment()  # raises if the manifest was edited after the run
    written: tuple[Path, ...] = ()
    if original.kind == "well_mixed":
        result = run_well_mixed(config)
    elif original.kind == "reactive_transport":
        if args.output:
            result, written = _run_in_space(config, Path(args.output))
        else:
            result = run_reactive_transport(config)
    elif original.kind != "ecosystem":
        result = run(config)
    elif args.output:
        # The replay's own outputs, streamed as it runs.
        result, _, written = _run_ecosystem_into(
            config, Path(args.output), _frame_every(config, None)
        )
    else:
        result = run_ecosystem(config, frame_every=None)  # only the final state is compared

    recorded = _comparable(original.outputs)
    fresh = _comparable(result.manifest.outputs)
    differences = [
        f"  {key}: recorded {recorded[key]!r}, replayed {value!r}"
        for key, value in fresh.items()
        if recorded.get(key) != value
    ]

    print(f"run_id      {original.run_id}")
    print(f"experiment  {original.experiment_id}")
    print(f"kind        {config.kind}")
    print(f"seed        {original.seed}")
    print(f"checksum    {original.config_sha256[:16]}... verified")
    if original.environment != result.manifest.environment:
        print("\nnote: the software environment differs from the original run")
        for key, was in sorted(original.environment.items()):
            now = result.manifest.environment.get(key)
            if was != now:
                print(f"  {key}: recorded {was}, now {now}")
    if differences:
        print("\nreplay DIFFERS from the recorded run:")
        print("\n".join(differences))
        if written:
            print(f"\nthe outputs written to {args.output} come from this differing replay")
        return 1
    print("\nreplay reproduced the recorded results exactly")
    if args.output and not isinstance(result, EcosystemResult | ReactiveTransportResult):
        written = _write_outputs(result, Path(args.output))
    for path in written:
        print(f"wrote {path}")
    return 0


def _command_ecosystem(args: argparse.Namespace) -> int:
    config = load_ecosystem_experiment(args.experiment)
    output_dir = (
        Path(args.output)
        if args.output
        else Path(args.experiment).parent / "runs" / config.experiment_id
    )
    frame_every = _frame_every(config, args.frame_interval_h)
    result, stored, (frames, viewer, manifest) = _run_ecosystem_into(
        config, output_dir, frame_every
    )
    print(f"run_id      {result.manifest.run_id}")
    print(f"experiment  {config.experiment_id}")
    print(f"steps       {result.final_state.step} over {result.final_state.time_h:g} h")
    print(f"frames      {stored} stored, one every {frame_every} step(s)")
    print(f"wrote       {frames}")
    print(f"wrote       {viewer}")
    print(f"wrote       {manifest}")
    print(f"\nreplay it with:  marse replay {manifest}")
    return 0


def _command_niche_scan(args: argparse.Namespace) -> int:
    scan = load_niche_scan(args.experiment)
    result = run_niche_scan(scan)
    output_dir = (
        Path(args.output)
        if args.output
        else Path(args.experiment).parent / "runs" / scan.experiment_id
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = result.write_json(output_dir / "niche-scan.json")
    csv_path = result.write_csv(output_dir / "niche-scan.csv")
    print(f"experiment  {scan.experiment_id}")
    print(f"rows        {len(result.rows)}")
    print(f"wrote       {json_path}")
    print(f"wrote       {csv_path}")
    return 0


def _command_ecosystem_batch(args: argparse.Namespace) -> int:
    raw = json.loads(Path(args.experiment).read_text(encoding="utf-8"))
    batch = ScenarioBatch.grid(raw["base"], raw.get("parameters", {}))
    catalog = run_batch(batch, workers=args.workers)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    catalog_path = catalog.write_json(output / "catalog.json")
    print(f"scenarios   {len(batch.scenarios)}")
    print(f"completed   {len(catalog.query(status='completed'))}")
    print(f"failed      {len(catalog.query(status='failed'))}")
    print(f"wrote       {catalog_path}")
    return 0


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
    if scales["ratio"] is not None:
        print(
            f"time scales diffusion across the box {scales['diffusion_h'] * 3600:.3g} s, "
            f"fastest process {scales['growth_h']:.3g} h (ratio {scales['ratio']:.2g})"
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
            brings = f"releases {what or 'nothing'} mmol"
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
    print(f"              the top voxel comes to it in about: {seconds}")
    if config.domain.mouth is not None:
        print("              it holds the mouth's saliva at saturation too")


def _report_plaque(config: ReactiveTransportConfig) -> None:
    from marse.biofilm.spreading import SPREADING_VERSION

    plaque = config.domain.plaque
    assert plaque is not None
    fills = ", ".join(f"{n} at {v:g} mol per m3" for n, v in plaque.packing_mol_per_m3.items())
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="marse",
        description="MARSE: Microbial Adaptability Resource Simulation Engine.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command")

    run_command = commands.add_parser("run", help="run a simulation from an experiment file")
    run_command.add_argument("experiment", help="path to a JSON experiment configuration")
    run_command.add_argument(
        "-o", "--output", help="directory for outputs (default: runs/<run_id>)"
    )
    run_command.set_defaults(handler=_command_run)

    replay_command = commands.add_parser("replay", help="re-run from a manifest and compare")
    replay_command.add_argument("manifest", help="path to a manifest.json written by 'marse run'")
    replay_command.add_argument("-o", "--output", help="directory for the replayed outputs")
    replay_command.set_defaults(handler=_command_replay)

    ecosystem_command = commands.add_parser(
        "ecosystem", help="run a 2D multi-species ecosystem and export a browser viewer"
    )
    ecosystem_command.add_argument("experiment", help="path to an ecosystem JSON configuration")
    ecosystem_command.add_argument("-o", "--output", help="directory for frames and viewer")
    ecosystem_command.add_argument(
        "--frame-interval-h",
        type=float,
        help=f"hours between stored frames (default: about {MAX_STORED_FRAMES} frames per run)",
    )
    ecosystem_command.set_defaults(handler=_command_ecosystem)

    niche_command = commands.add_parser(
        "niche-scan", help="scan environmental conditions against species capabilities"
    )
    niche_command.add_argument("experiment", help="path to a niche scan JSON configuration")
    niche_command.add_argument("-o", "--output", help="directory for scan outputs")
    niche_command.set_defaults(handler=_command_niche_scan)

    batch_command = commands.add_parser(
        "ecosystem-batch", help="run a reproducible ecosystem parameter ensemble"
    )
    batch_command.add_argument("experiment", help="JSON file with base and parameters")
    batch_command.add_argument("-o", "--output", required=True, help="directory for catalog output")
    batch_command.add_argument("--workers", type=int, default=1, help="parallel worker processes")
    batch_command.set_defaults(handler=_command_ecosystem_batch)

    check_command = commands.add_parser(
        "check", help="check a reaction network and print its balanced processes"
    )
    check_command.add_argument("network", help="path to a network JSON file (schema version 2)")
    check_command.set_defaults(handler=_command_check)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return 0
    errors = np.geterr()
    try:
        return int(args.handler(args))
    except ConfigError as error:
        print(f"marse: invalid experiment: {error}")
        return 2
    except NicheError as error:
        print(f"marse: invalid niche scan: {error}")
        return 2
    except (OSError, ValueError, KeyError) as error:
        print(f"marse: {type(error).__name__}: {error}")
        return 2
    finally:
        np.seterr(**errors)  # NumPy's error handling as the caller had it
