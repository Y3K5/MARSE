"""Solving the oxygen profile through a biofilm, and where it runs out.

Where examples/oxygen_penetration.py estimated the penetration depth from a
closed-form expression, this one solves the reaction-diffusion equation
numerically and compares the two (docs/theory.md, section 5; validation case
V3):

    D d2C/dz2 = k_max C / (K + C),   C(0) = bulk,   dC/dz = 0 at the base.

Run with:  python examples/biofilm_profile.py

Uptake rates and biomass densities are illustrative (confidence C in
docs/parameters.md).
"""

from __future__ import annotations

import numpy as np

from marse.spatial.diffusion import solve_steady_state
from marse.spatial.domain import Grid1D
from marse.spatial.solutes import (
    OXYGEN_MOLAR_MASS,
    oxygen_diffusivity_m2_per_s,
    oxygen_saturation_mg_per_l,
)
from marse.validation.analytical import zero_order_penetration_depth

RELATIVE_DIFFUSIVITY = 0.43  # small nonpolar solutes in biofilms (Stewart 1998)
THICKNESS_UM = 200.0
CELLS = 800
BIOMASS_G_PER_L = 25.0
SPECIFIC_UPTAKE = 24.0  # mmol O2 per g dry weight per h
HALF_SATURATION_MM = 1.0e-3  # oxygen affinity, illustrative


def bar(value: float, reference: float, width: int = 32) -> str:
    filled = round(width * max(value, 0.0) / reference) if reference > 0 else 0
    return "#" * filled + "." * (width - filled)


def main() -> None:
    print("MARSE example: oxygen profile through a biofilm")
    print("=" * 64)

    surface_mm = float(oxygen_saturation_mg_per_l(37.0)) / OXYGEN_MOLAR_MASS
    diffusivity = RELATIVE_DIFFUSIVITY * float(oxygen_diffusivity_m2_per_s(37.0)) * 1e12
    max_uptake = BIOMASS_G_PER_L * SPECIFIC_UPTAKE / 3600.0  # mM/s

    grid = Grid1D(thickness=THICKNESS_UM, cells=CELLS)
    print(f"\nbiofilm       {THICKNESS_UM:.0f} um thick, {CELLS} nodes ({grid.dx:.2f} um apart)")
    print(f"surface O2    {surface_mm:.4f} mM (air-saturated water at 37 C)")
    print(
        f"diffusivity   {diffusivity:.0f} um2/s ({RELATIVE_DIFFUSIVITY:.2f} of the value in water)"
    )
    print(f"max uptake    {max_uptake:.4f} mM/s at {BIOMASS_G_PER_L:.0f} g/L")

    profile = solve_steady_state(
        grid,
        diffusivity=diffusivity,
        surface=surface_mm,
        max_uptake=max_uptake,
        half_saturation=HALF_SATURATION_MM,
    )
    print(f"\nsolved in {profile.iterations} Newton iterations")

    print(f"\n  {'z (um)':>7}  {'C (mM)':>9}  {'% surface':>9}  profile")
    for z in np.linspace(0.0, THICKNESS_UM, 17)[1:]:
        i = int(np.argmin(np.abs(profile.depths - z)))
        c = profile.concentration[i]
        print(
            f"  {profile.depths[i]:>7.1f}  {c:>9.5f}  {100 * c / surface_mm:>8.1f}%  "
            f"{bar(c, surface_mm)}"
        )

    measured = profile.penetration_depth()
    analytical = zero_order_penetration_depth(diffusivity, surface_mm, max_uptake)
    print(f"\npenetration depth (1% of surface)   {measured:.1f} um   [solved]")
    print(f"zero-order estimate  sqrt(2 D C0/k)  {analytical:.1f} um   [closed form]")
    print("  the closed form marks where C reaches zero; at a 1% threshold the")
    print(f"  parabola gives 0.9 x that, so {0.9 * analytical:.1f} um is the like-for-like value.")
    print(f"\nanoxic fraction of the biofilm       {profile.anoxic_fraction:.1%}")

    flux = profile.surface_flux(diffusivity)
    uptake = profile.total_uptake(max_uptake, HALF_SATURATION_MM)
    print("\nconservation check")
    print(f"  oxygen entering the surface  {flux:.6f} mM um/s")
    print(f"  oxygen consumed inside       {uptake:.6f} mM um/s")
    print(f"  relative difference          {abs(flux - uptake) / uptake:.2e}")

    print("\nMost of this biofilm is anoxic. That is the result of theory.md")
    print("section 5.1: penetration scales as the square root of supply over")
    print("demand, so an active biofilm thicker than a few tens of micrometres")
    print("is stratified no matter how the parameters are adjusted.")


if __name__ == "__main__":
    main()
