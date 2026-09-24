"""Batch growth of a single species under substrate and temperature limitation.

Demonstrates, with no simulation kernel yet, that the pieces MARSE is built
from already agree with each other:

  * Monod growth coupled to substrate uptake (docs/theory.md, sections 1.3, 3.2)
  * the exact integrated batch solution as an independent reference (section 3.5)
  * temperature scaling through the cardinal model (section 2.2)
  * conservation of mass, checked rather than assumed (section 3.4)

Run with:  python examples/batch_growth.py

Parameters are illustrative (confidence C in docs/parameters.md); they are here
to make the example runnable, not to describe a particular organism.
"""

from __future__ import annotations

import numpy as np

from marse.microbes.cardinal import cardinal_temperature
from marse.microbes.growth import doubling_time, monod, substrate_uptake_rate
from marse.validation.analytical import batch_final_biomass, monod_batch_time

# Illustrative parameters; units are hours, mM and g/L throughout.
MU_OPT = 1.0  # optimal specific growth rate, 1/h
K_S = 0.5  # half-saturation constant, mM
# 0.45 g biomass per g glucose (docs/parameters.md) at 180.16 g/mol:
YIELD = 0.45 * 180.16 / 1000.0  # biomass yield, g dry weight per mmol
MAINTENANCE = 0.0  # maintenance coefficient, mmol per g per h
S0 = 22.2  # initial substrate, mM (0.4 % w/v glucose)
X0 = 0.01  # initial biomass, g/L
T_CARDINAL = (6.0, 40.0, 47.0)  # T_min, T_opt, T_max in degrees Celsius


def simulate(mu_max: float, t_end: float, steps: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Integrate the coupled batch system with fourth-order Runge-Kutta."""

    def rhs(state: np.ndarray) -> np.ndarray:
        s, x = state
        mu = monod(s, mu_max, K_S)
        uptake = substrate_uptake_rate(mu, x, YIELD, MAINTENANCE)
        return np.array([-uptake, mu * x])

    h = t_end / steps
    state = np.array([S0, X0])
    trajectory = [state]
    for _ in range(steps):
        k1 = rhs(state)
        k2 = rhs(state + h / 2 * k1)
        k3 = rhs(state + h / 2 * k2)
        k4 = rhs(state + h * k3)
        state = state + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        trajectory.append(state)
    history = np.array(trajectory)
    return np.linspace(0.0, t_end, steps + 1), history[:, 0], history[:, 1]


def main() -> None:
    print("MARSE example: batch growth under substrate limitation")
    print("=" * 62)

    # --- how temperature scales the maximum growth rate ---------------------
    print("\nTemperature scaling (cardinal model, T_min/T_opt/T_max = 6/40/47 C)")
    print(f"  {'T (C)':>7}  {'gamma_T':>8}  {'mu_max (1/h)':>13}  {'doubling (min)':>15}")
    for temperature in (10.0, 20.0, 30.0, 37.0, 40.0, 45.0):
        gamma = float(cardinal_temperature(temperature, *T_CARDINAL))
        mu_max = MU_OPT * gamma
        doubling = f"{doubling_time(mu_max) * 60:.1f}" if mu_max > 0 else "no growth"
        print(f"  {temperature:>7.1f}  {gamma:>8.3f}  {mu_max:>13.3f}  {doubling:>15}")

    # --- a batch culture at 37 C -------------------------------------------
    mu_max = MU_OPT * float(cardinal_temperature(37.0, *T_CARDINAL))
    t, substrate, biomass = simulate(mu_max, t_end=25.0, steps=25_000)

    print(f"\nBatch culture at 37 C (mu_max = {mu_max:.3f} 1/h)")
    print(f"  {'t (h)':>6}  {'S (mM)':>9}  {'X (g/L)':>9}  {'mu (1/h)':>9}")
    for hour in (0.0, 2.0, 4.0, 6.0, 7.0, 8.0, 9.0, 12.0, 25.0):
        i = int(hour / t[-1] * (len(t) - 1))
        mu = float(monod(substrate[i], mu_max, K_S))
        print(f"  {t[i]:>6.1f}  {substrate[i]:>9.4f}  {biomass[i]:>9.4f}  {mu:>9.4f}")

    # --- checks against the analytical reference ----------------------------
    print("\nChecks")
    expected_final = batch_final_biomass(S0, X0, YIELD)
    print(f"  final biomass, simulated : {biomass[-1]:.6f} g/L")
    print(f"  final biomass, predicted : {expected_final:.6f} g/L  (X0 + Y*S0)")

    conserved = biomass + YIELD * substrate
    drift = float(np.max(np.abs(conserved - conserved[0])) / conserved[0])
    print(f"  mass-balance drift       : {drift:.2e}  (X + Y*S must stay constant)")

    # The exact integrated solution: time to reach a given substrate level.
    print(f"\n  {'S (mM)':>8}  {'t simulated':>12}  {'t analytical':>13}  {'rel. error':>11}")
    for target in (15.0, 10.0, 5.0, 1.0, 0.5):
        i = int(np.argmin(np.abs(substrate - target)))
        exact = float(monod_batch_time(substrate[i], S0, X0, mu_max, K_S, YIELD))
        error = abs(t[i] - exact) / max(exact, 1e-12)
        print(f"  {substrate[i]:>8.3f}  {t[i]:>12.5f}  {exact:>13.5f}  {error:>11.2e}")

    print("\nThe simulated trajectory reproduces the closed-form solution, and")
    print("mass is conserved to machine precision. See docs/theory.md sections")
    print("1.3, 3.2-3.5 for the equations and their assumptions.")


if __name__ == "__main__":
    main()
