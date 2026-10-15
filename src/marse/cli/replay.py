"""``marse replay``: re-run a manifest and compare it with the recorded result."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from marse.cli.present import warn_unverified
from marse.core.provenance import Manifest
from marse.engines import registry

_BIOFILM_KEYS = ("penetration_depth_um", "active_zone_um", "base_concentration_mm")


def command_replay(args: argparse.Namespace) -> int:
    original = Manifest.read(args.manifest)
    engine = registry.for_kind(original.kind)
    if engine.unverified:
        warn_unverified()
    output = Path(args.output) if args.output else None
    completed = engine.replay(original, output if engine.writes_during_replay else None)
    result = completed.result
    written = completed.written
    fresh_manifest = result.manifest
    if fresh_manifest is None:
        raise ValueError("replay produced no manifest to compare")

    recorded = _comparable(original.outputs)
    fresh = _comparable(fresh_manifest.outputs)
    differences = [
        f"  {key}: recorded {recorded[key]!r}, replayed {value!r}"
        for key, value in fresh.items()
        if recorded.get(key) != value
    ]

    print(f"run_id      {original.run_id}")
    print(f"experiment  {original.experiment_id}")
    print(f"kind        {result.config.kind}")
    print(f"seed        {original.seed}")
    print(f"checksum    {original.config_sha256[:16]}... verified")
    if original.environment != fresh_manifest.environment:
        print("\nnote: the software environment differs from the original run")
        for key, was in sorted(original.environment.items()):
            now = fresh_manifest.environment.get(key)
            if was != now:
                print(f"  {key}: recorded {was}, now {now}")
    if differences:
        print("\nreplay DIFFERS from the recorded run:")
        print("\n".join(differences))
        if written:
            print(f"\nthe outputs written to {args.output} come from this differing replay")
        return 1
    print("\nreplay reproduced the recorded results exactly")
    if output is not None and not engine.writes_during_replay:
        written = engine.write_outputs(result, output)
    for path in written:
        print(f"wrote {path}")
    return 0


def _comparable(outputs: dict[str, Any]) -> dict[str, float | str]:
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
        if "spreading" in outputs:
            final = outputs["spreading"]["final"]
            values.update({f"structure[{k}]": float(v) for k, v in final.items()})
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
