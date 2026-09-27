"""How many early colonizers stick to four dental surfaces, and how fast.

Runs the dental scene of examples/environments/dental/dental_surfaces.json as
one thin column per material, enamel, titanium, zirconia and acrylic (PMMA),
all under a salivary pellicle and the same salivary film. Streptococcus oralis
and S. sanguinis bind from the saliva, lock, and grow on salivary glucose. The
table gives bound cells per cm² and the fraction of the surface under cells.

The only contrast between the materials that data supports is titanium against
zirconia (Scarano et al. 2004; Oda et al. 2020); the other parameters are
illustrative (docs/parameters.md, "Surfaces and adhesion").

Run with: python examples/surface_adhesion.py
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

from marse.core.reactive_transport import run
from marse.schemas.experiment import experiment_from_dict

SCENE = Path(__file__).resolve().parent / "environments" / "dental" / "dental_surfaces.json"
HOURS = (1, 2, 6, 12, 24)


def column(raw: dict, material: str) -> dict:
    """The scene as one column of voxels over a single material."""
    raw = copy.deepcopy(raw)
    domain = raw["domain"]
    domain["voxels"] = domain["voxels"][-1:]
    domain["substratum"]["patches"] = [{"material": material, "region_um": []}]
    domain["adhesion"] = [a for a in domain["adhesion"] if a["material"] == material]
    return raw


def main() -> None:
    raw = json.loads(SCENE.read_text("utf-8"))
    materials = [p["material"] for p in raw["domain"]["substratum"]["patches"]]
    species = [s["attached"] for s in raw["domain"]["suspension"]]
    print("bound cells per cm2 (all species) and the fraction covered, by hour")
    print(f"{'material':<10}" + "".join(f"{h:>12} h" for h in HOURS))
    for material in materials:
        result = run(experiment_from_dict(column(raw, material)))
        cells, covered = [], []
        for hour in HOURS:
            i = list(result.times_h).index(float(hour))
            cells.append(sum(result.surface[f"{material}_{s}_cells_per_cm2"][i] for s in species))
            covered.append(result.surface[f"{material}_covered"][i])
        print(f"{material:<10}" + "".join(f"{n:>14.3g}" for n in cells))
        print(f"{'':<10}" + "".join(f"{c:>14.1%}" for c in covered))


if __name__ == "__main__":
    main()
