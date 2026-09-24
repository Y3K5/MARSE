"""How far oxygen reaches into a biofilm, and why that makes biofilms stratified.

Demonstrates (docs/theory.md, sections 4.3-4.5 and 5):

  * temperature-dependent oxygen solubility and diffusivity
  * the reduction of diffusivity inside a biofilm matrix
  * the zero-order penetration depth, delta = sqrt(2 D_e C0 / k0)
  * the separation of time scales that justifies a quasi-steady solute field
  * the explicit-solver stability limit that makes that approximation necessary

Run with:  python examples/oxygen_penetration.py

Uptake rates and biomass densities here are illustrative (confidence C in
docs/parameters.md).
"""

from __future__ import annotations

import numpy as np

from marse.spatial.solutes import (
    OXYGEN_MOLAR_MASS,
    oxygen_diffusivity_m2_per_s,
    oxygen_saturation_mg_per_l,
)
from marse.validation.analytical import zero_order_penetration_depth

RELATIVE_DIFFUSIVITY = 0.43  # small nonpolar solutes in biofilms (Stewart 1998)
BIOMASS_DENSITY = 25.0  # g dry weight per L, illustrative for a dense biofilm
SPECIFIC_UPTAKE = 24.0  # mmol O2 per g dry weight per h, illustrative
GRID_SPACING_UM = 5.0  # a typical spatial resolution


def to_micrometre_units(diffusivity_m2_per_s: float) -> float:
    """Convert m^2/s to um^2/s, the internal length scale for a biofilm."""
    return diffusivity_m2_per_s * 1e12


def main() -> None:
    print("MARSE example: oxygen penetration into a biofilm")
    print("=" * 62)

    # --- oxygen availability and mobility versus temperature ----------------
    print("\nOxygen in water at equilibrium with air (1 atm)")
    print(f"  {'T (C)':>6}  {'C* (mg/L)':>10}  {'C* (mM)':>9}  {'D_water (um2/s)':>16}")
    for temperature in (20.0, 25.0, 30.0, 37.0):
        saturation = float(oxygen_saturation_mg_per_l(temperature))
        diffusivity = to_micrometre_units(float(oxygen_diffusivity_m2_per_s(temperature)))
        print(
            f"  {temperature:>6.1f}  {saturation:>10.2f}  "
            f"{saturation / OXYGEN_MOLAR_MASS:>9.3f}  {diffusivity:>16.0f}"
        )
    print("  note: warming water holds LESS oxygen but diffuses it FASTER;")
    print("        the two effects oppose each other (theory.md 4.4-4.5).")

    # --- conditions inside the biofilm at 37 C ------------------------------
    c0 = float(oxygen_saturation_mg_per_l(37.0)) / OXYGEN_MOLAR_MASS  # mM
    d_water = to_micrometre_units(float(oxygen_diffusivity_m2_per_s(37.0)))
    d_effective = RELATIVE_DIFFUSIVITY * d_water
    uptake = BIOMASS_DENSITY * SPECIFIC_UPTAKE / 3600.0  # mM/s

    print("\nInside the biofilm at 37 C")
    print(f"  surface oxygen concentration  C0  = {c0:.4f} mM")
    print(f"  diffusivity in water          D   = {d_water:.0f} um2/s")
    print(f"  effective diffusivity  ({RELATIVE_DIFFUSIVITY:.2f} D)  D_e = {d_effective:.0f} um2/s")
    print(f"  volumetric uptake             k0  = {uptake:.4f} mM/s")

    depth = float(zero_order_penetration_depth(d_effective, c0, uptake))
    print(f"\n  penetration depth  delta = sqrt(2 D_e C0 / k0) = {depth:.1f} um")

    # --- the resulting profile ----------------------------------------------
    print("\n  Steady-state profile, C(z) = C0 (1 - z/delta)^2:")
    print(f"  {'z (um)':>8}  {'C (mM)':>9}  {'% of C0':>8}")
    for z in np.linspace(0.0, depth * 1.4, 8):
        concentration = c0 * (1 - z / depth) ** 2 if z < depth else 0.0
        state = "" if z < depth else "   <- anoxic"
        print(f"  {z:>8.1f}  {concentration:>9.4f}  {concentration / c0 * 100:>7.1f}%{state}")

    # --- why a thick biofilm is necessarily stratified ----------------------
    print("\nPenetration depth versus demand (delta scales as 1/sqrt(k0))")
    print(f"  {'biomass (g/L)':>14}  {'k0 (mM/s)':>10}  {'delta (um)':>11}")
    for density in (5.0, 10.0, 25.0, 50.0, 100.0):
        k = density * SPECIFIC_UPTAKE / 3600.0
        print(
            f"  {density:>14.0f}  {k:>10.4f}  "
            f"{float(zero_order_penetration_depth(d_effective, c0, k)):>11.1f}"
        )
    print("  quadrupling the demand only halves the depth, so an active biofilm")
    print("  thicker than a few tens of um is always stratified (theory.md 5.1).")

    # --- time scales ---------------------------------------------------------
    thickness = 100.0
    tau_diffusion = thickness**2 / d_effective
    print(f"\nTime scales for a {thickness:.0f} um biofilm")
    print(f"  solute equilibration  tau_D = L^2/D_e = {tau_diffusion:.1f} s")
    for doubling_minutes in (20.0, 60.0):
        tau_growth = doubling_minutes * 60.0
        print(
            f"  biomass doubling ({doubling_minutes:.0f} min) = {tau_growth:.0f} s"
            f"   ->  ratio {tau_growth / tau_diffusion:.0f}x slower"
        )
    print("  solutes equilibrate orders of magnitude faster than biomass changes,")
    print("  which is what justifies a quasi-steady solute field (theory.md 5.3).")

    # --- and why the field must be solved, not stepped ----------------------
    max_step = GRID_SPACING_UM**2 / (4 * d_water)
    print(f"\nExplicit-solver stability limit (2D, {GRID_SPACING_UM:.0f} um grid)")
    print(f"  dt <= h^2 / (2 d D) = {max_step * 1e3:.2f} ms")
    print(f"  stepping one hour would take {3600 / max_step:,.0f} steps for the")
    print("  solute field alone, so MARSE solves it to steady state instead.")


if __name__ == "__main__":
    main()
