"""Command-line entry point: ``marse``.

Two commands matter at this stage:

    marse run experiment.json        run a simulation and write its manifest
    marse replay manifest.json       re-run from a manifest and compare

``replay`` is the reproducibility claim made executable: it rebuilds the
configuration from the manifest, verifies the checksum, runs again, and reports
whether the results match.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from marse import __version__
from marse.core.config import ConfigError, load_experiment
from marse.core.provenance import Manifest
from marse.core.simulation import BiofilmProfileResult, SimulationResult, run

Result = SimulationResult | BiofilmProfileResult


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


def _comparable(outputs: dict) -> dict[str, float]:
    """The numbers a replay must reproduce, flattened so the two kinds compare alike."""
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


def _command_run(args: argparse.Namespace) -> int:
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


def _command_replay(args: argparse.Namespace) -> int:
    original = Manifest.read(args.manifest)
    config = original.experiment()  # raises if the manifest was edited after the run
    result = run(config)

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
        return 1
    print("\nreplay reproduced the recorded results exactly")
    if args.output:
        trajectory, manifest = _write_outputs(result, Path(args.output))
        print(f"wrote {trajectory}")
        print(f"wrote {manifest}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="marse",
        description="MARSE: Microbial Adaptability Resource Engine.",
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
    except (OSError, ValueError, KeyError) as error:
        print(f"marse: {type(error).__name__}: {error}")
        return 2
    finally:
        np.seterr(all="warn")
