"""Run controls and environmental variants for the periodontal example.

This is a compact, provenance-aware screening workflow, not a parameter
fitting or clinical prediction tool. Moisture is represented as an explicit
dimensionless condition and also scales available carbon/oxygen in the
illustrative setup; replace this mapping with measured medium data before
biological interpretation.

Run:
    python examples/periodontal_variant_analysis.py --output /tmp/periodontal-analysis
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

from marse.ecosystem.model import ecosystem_from_dict, run

BASE = Path(__file__).parent / "experiments" / "periodontal_pathogen_biofilm.json"
PROVENANCE = Path(__file__).parent / "experiments" / "periodontal_pathogen_provenance.json"
TEMPERATURES = (30.0, 37.0, 40.0)
MOISTURES = ("reduced", "reference")
REPLICATES = (0, 1, 2)


def _checksum(value: Any) -> str:
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _variant(
    base: dict[str, Any],
    species_names: list[str],
    temperature: float,
    moisture: str,
    seed: int,
    competition_mode: str = "reference",
    resource_perturbation: str = "reference",
) -> dict[str, Any]:
    config = deepcopy(base)
    label = species_names[0] if len(species_names) == 1 else "mixed"
    config["experiment_id"] = f"periodontal-{label}-{temperature:g}C-{moisture}-r{seed}"
    config["seed"] = seed
    config["duration_h"] = 6.0
    config["species"] = [item for item in config["species"] if item["name"] in species_names]
    count = len(config["species"])
    config["competition_coefficients"] = [
        [
            0.45 if row == column else (0.8 if competition_mode == "reference" else 0.0)
            for column in range(count)
        ]
        for row in range(count)
    ]
    config["conditions"][0]["initial"] = temperature
    moisture_factor = 0.65 if moisture == "reduced" else 1.0
    config["conditions"][1]["initial"] = moisture_factor
    for nutrient in config["nutrients"]:
        if nutrient["name"] in {"carbon", "oxygen"}:
            nutrient["initial"] *= moisture_factor
            if nutrient.get("boundary_value") is not None:
                nutrient["boundary_value"] *= moisture_factor
        if nutrient["name"] == resource_perturbation:
            nutrient["initial"] *= 0.5
            if nutrient.get("boundary_value") is not None:
                nutrient["boundary_value"] *= 0.5
    return config


def _record(config: dict[str, Any]) -> dict[str, Any]:
    result = run(ecosystem_from_dict(config))
    final = result.final_state
    record: dict[str, Any] = {
        "scenario_id": _checksum(config)[:16],
        "experiment_id": config["experiment_id"],
        "temperature_c": config["conditions"][0]["initial"],
        "moisture": config["conditions"][1]["initial"],
        "species_set": "mixed" if len(config["species"]) > 1 else config["species"][0]["name"],
        "seed": config["seed"],
    }
    for index, species in enumerate(config["species"]):
        record[f"biomass:{species['name']}"] = float(final.biomass[index].sum())
        record[f"occupied:{species['name']}"] = int(
            (final.biomass[index] > 0.01 * config["carrying_capacity"]).sum()
        )
    for index, nutrient in enumerate(config["nutrients"]):
        record[f"nutrient_mean:{nutrient['name']}"] = float(final.nutrients[index].mean())
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = json.loads(BASE.read_text(encoding="utf-8"))
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    species = [item["name"] for item in base["species"]]
    records = []
    for temperature in TEMPERATURES:
        for moisture in MOISTURES:
            for replicate in REPLICATES:
                for selected in (species, *([name] for name in species)):
                    records.append(
                        _record(
                            _variant(
                                base,
                                list(selected),
                                temperature,
                                moisture,
                                replicate + 1000,
                            )
                        )
                    )
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = sorted({key for record in records for key in record})
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    manifest = {
        "base_experiment": str(BASE.relative_to(BASE.parents[1])),
        "base_sha256": _checksum(base),
        "provenance_file": str(PROVENANCE.relative_to(PROVENANCE.parents[1])),
        "provenance_sha256": _checksum(provenance),
        "provenance_status": provenance["status"],
        "parameter_group_ids": [group["id"] for group in provenance["parameter_groups"]],
        "calibration_status": sorted(
            {group["calibration_status"] for group in provenance["parameter_groups"]}
        ),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "temperatures_c": TEMPERATURES,
        "moisture_levels": MOISTURES,
        "replicates": REPLICATES,
        "scenario_count": len(records),
        "notes": "Illustrative moisture mapping; not a measured periodontal moisture model.",
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"wrote {len(records)} scenario summaries to {args.output}")


if __name__ == "__main__":
    main()
