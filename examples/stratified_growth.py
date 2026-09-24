"""Why a thick biofilm grows no faster than a thin one.

Couples the solute field to the biomass that consumes it (docs/theory.md,
sections 5 and 6): each depth sees its own oxygen concentration and therefore
grows at its own rate. The consequence is stratification — past the
penetration depth, added biomass adds no growth, so the biofilm-averaged rate
falls as the film thickens even though every cell is identical.

Run with:  python examples/stratified_growth.py

Parameters are illustrative (confidence C in docs/parameters.md). Units are
hours and micrometres throughout, which is why the diffusivity is fetched in
um^2/h rather than the um^2/s that reference tables quote.
"""

from __future__ import annotations

import numpy as np

from marse.biofilm.biomass import BiofilmState, Population, solve_growth_profile
from marse.spatial.domain import Grid1D
from marse.spatial.solutes import (
    OXYGEN_MOLAR_MASS,
    oxygen_diffusivity_um2_per_h,
    oxygen_saturation_mg_per_l,
)

RELATIVE_DIFFUSIVITY = 0.43  # small nonpolar solutes in biofilms (Stewart 1998)
DENSITY_G_PER_L = 25.0
AEROBE = Population("aerobe", mu_max=0.3, half_saturation=1.0e-3, yield_per_substrate=0.08)


def bar(value: float, reference: float, width: int = 30) -> str:
    filled = round(width * max(value, 0.0) / reference) if reference > 0 else 0
    return "#" * filled + "." * (width - filled)


def main() -> None:
    print("MARSE example: growth stratified by an oxygen gradient")
    print("=" * 68)

    surface = float(oxygen_saturation_mg_per_l(37.0)) / OXYGEN_MOLAR_MASS
    diffusivity = RELATIVE_DIFFUSIVITY * float(oxygen_diffusivity_um2_per_h(37.0))
    print(f"\nsurface oxygen  {surface:.4f} mM        biomass  {DENSITY_G_PER_L:.0f} g/L, uniform")
    print(
        f"diffusivity     {diffusivity:.3e} um2/h  ({RELATIVE_DIFFUSIVITY} of the value in water)"
    )
    print(f"mu_max          {AEROBE.mu_max:.2f} /h            K  {AEROBE.half_saturation:g} mM")

    # --- the profile through one thick biofilm ------------------------------
    thickness = 400.0
    grid = Grid1D(thickness, 1200)
    state = BiofilmState.uniform(grid, (AEROBE,), [DENSITY_G_PER_L])
    result = solve_growth_profile(state, diffusivity=diffusivity, surface=surface)
    peak = float(np.max(result.growth_rate[0]))

    print(f"\nA {thickness:.0f} um biofilm, solved in {result.solute.iterations} iterations:")
    print(f"\n  {'z (um)':>7}  {'O2 (mM)':>9}  {'mu (1/h)':>9}  {'% of peak':>9}  growth")
    for z in np.linspace(0.0, thickness, 17)[1:]:
        i = int(np.argmin(np.abs(result.depths - z)))
        mu = float(result.growth_rate[0][i])
        print(
            f"  {result.depths[i]:>7.1f}  {result.solute.concentration[i]:>9.5f}  "
            f"{mu:>9.5f}  {100 * mu / peak:>8.1f}%  {bar(mu, peak)}"
        )

    print(f"\n  oxygen penetrates      {result.solute.penetration_depth():>6.1f} um")
    print(f"  active zone (>50% mu)  {result.active_zone():>6.1f} um")
    print(f"  of the biofilm         {result.active_zone() / thickness:>6.1%}")

    # --- thicker films do not grow faster -----------------------------------
    print("\nThe same biofilm at different thicknesses:")
    print(f"\n  {'thickness':>10}  {'active zone':>12}  {'mean mu':>9}  {'vs surface':>10}")
    reference = None
    for total in (25.0, 50.0, 100.0, 200.0, 400.0, 800.0):
        grid = Grid1D(total, max(100, int(total * 3)))
        state = BiofilmState.uniform(grid, (AEROBE,), [DENSITY_G_PER_L])
        run = solve_growth_profile(state, diffusivity=diffusivity, surface=surface)
        mean = float(run.mean_growth_rate()[0])
        surface_rate = float(run.growth_rate[0][0])
        reference = reference or mean
        print(
            f"  {total:>9.0f}um  {run.active_zone():>11.1f}um  {mean:>9.5f}  "
            f"{mean / surface_rate:>9.1%}"
        )

    print("\nThe active zone stops deepening once the film passes the penetration")
    print("depth, so from there on the mean growth rate roughly halves whenever")
    print("the thickness doubles. Every cell is identical; the gradient alone")
    print("decides which of them grow. This is the mechanism behind the narrow")
    print("bands of protein synthesis reported in real biofilms, and behind the")
    print("tolerance of their interiors to agents that only affect active cells.")


if __name__ == "__main__":
    main()
