"""Run periodontal single-, pairwise-, and mixed-community controls.

The output is deliberately compact: no trajectory frames are retained. This
study separates community composition, competition assumptions, and resource
limitation before any organism-specific calibration is attempted.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

try:
    from examples.periodontal_variant_analysis import (
        BASE,
        PROVENANCE,
        REPLICATES,
        _checksum,
        _variant,
    )
except ModuleNotFoundError:
    from periodontal_variant_analysis import BASE, PROVENANCE, REPLICATES, _checksum, _variant

from marse.ecosystem.model import ecosystem_from_dict, run

SPECIES = ("T. denticola", "T. forsythia", "P. gingivalis")
COMMUNITIES = (
    ("T. denticola",),
    ("T. forsythia",),
    ("P. gingivalis",),
    ("T. denticola", "T. forsythia"),
    ("T. denticola", "P. gingivalis"),
    ("T. forsythia", "P. gingivalis"),
    SPECIES,
)
RESOURCE_PERTURBATIONS = (
    "reference",
    "oxygen",
    "heme",
    "peptides",
    "succinate_ablation",
    "acetate_ablation",
)
COMPETITION_MODES = ("reference", "none")


def _record(config: dict, competition_mode: str, resource_perturbation: str) -> dict[str, object]:
    result = run(ecosystem_from_dict(config))
    final = result.final_state
    occupied = final.biomass > 0.01 * config["carrying_capacity"]
    record: dict[str, object] = {
        "scenario_id": _checksum(config)[:16],
        "community": "+".join(item["name"] for item in config["species"]),
        "competition_mode": competition_mode,
        "resource_perturbation": resource_perturbation,
        "replicate": config["seed"],
        "shared_occupied_cells": int((occupied.sum(axis=0) >= 2).sum()),
    }
    for index, species in enumerate(config["species"]):
        record[f"biomass:{species['name']}"] = float(final.biomass[index].sum())
        record[f"occupied:{species['name']}"] = int(occupied[index].sum())
    for index, additive in enumerate(config.get("additives", [])):
        record[f"metabolite_mean:{additive['name']}"] = float(final.additives[index].mean())
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = json.loads(BASE.read_text(encoding="utf-8"))
    records = []
    for competition_mode in COMPETITION_MODES:
        for resource in RESOURCE_PERTURBATIONS:
            for replicate in REPLICATES:
                for community in COMMUNITIES:
                    config = _variant(
                        base,
                        list(community),
                        37.0,
                        "reference",
                        2000 + replicate,
                        competition_mode,
                        resource,
                    )
                    records.append(_record(config, competition_mode, resource))
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = sorted({key for record in records for key in record})
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    manifest = {
        "base_experiment": str(BASE.relative_to(BASE.parents[1])),
        "base_sha256": _checksum(base),
        "provenance_sha256": _checksum(provenance),
        "study": "periodontal-community-controls",
        "temperature_c": 37.0,
        "moisture": "reference",
        "communities": ["+".join(community) for community in COMMUNITIES],
        "competition_modes": COMPETITION_MODES,
        "resource_perturbations": RESOURCE_PERTURBATIONS,
        "replicates": REPLICATES,
        "scenario_count": len(records),
        "interpretation": "Exploratory controls; parameters remain uncalibrated placeholders.",
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {len(records)} control summaries to {args.output}")


if __name__ == "__main__":
    main()
