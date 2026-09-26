"""Command-line entry point: ``marse``.

Four commands matter at this stage:

    marse run experiment.json        run a batch or biofilm simulation, or a
                                     version 2 reaction network in a closed box
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
from marse.core.config import ConfigError, load_experiment
from marse.core.provenance import Manifest
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
from marse.schemas.experiment import load_experiment as load_network_experiment
from marse.schemas.network import RUN_FIELDS, SCHEMA_VERSION, read_document

Result = SimulationResult | BiofilmProfileResult | WellMixedResult


def _write_outputs(result: Result, output_dir: Path) -> tuple[Path, Path]:
    """Write the run's data file and its manifest. The data file differs by kind."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if isinstance(result, BiofilmProfileResult):
        data = result.write_profile(output_dir / "profile.csv")
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
    worst = max(q["largest_relative_residual"] for q in outputs["balance"].values())
    print(f"balance     carbon, nitrogen and electrons conserved to {worst:.1e} of their totals")


def _comparable(outputs: dict) -> dict[str, float | str]:
    """The values a replay must reproduce, flattened so every kind compares alike."""
    if "final_mol_per_m3" in outputs:  # a reaction network in a closed box
        values: dict[str, float | str] = {
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


def _run_network(args: argparse.Namespace) -> int:
    config = load_network_experiment(args.experiment)
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
    if args.output and not isinstance(result, EcosystemResult):
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
    print(f"  {'name':<24} {'phase':<12} {'formula':<16} {'C':>6} {'N':>6} {'e-':>6} {'g/mol':>9}")
    for c in network.components:
        carbon, nitrogen, electrons = (f"{float(q):g}" for q in c.formula.composition)
        print(
            f"  {c.name:<24} {c.phase:<12} {c.formula.label():<16} {carbon:>6} {nitrogen:>6} "
            f"{electrons:>6} {c.formula.molar_mass_g_per_mol:>9.3f}"
        )
    print(
        f"\nprocesses   {len(network.processes)}, each conserving carbon, nitrogen and "
        "electrons exactly"
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
            print(
                f"\nrunnable    {config.duration_h:g} h in steps of at most "
                f"{config.timestep_h:g} h: marse run {path.name}"
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
        np.seterr(all="warn")
