# Chemotaxis

Species can optionally move up or down a declared nutrient, condition, or
additive gradient using `chemotaxis_field` and
`chemotaxis_sensitivity`. The default sensitivity is zero, so existing
experiments are unchanged.

The current implementation is a conservative finite-difference transport
term with no-flux outer edges. Positive sensitivity moves biomass toward
increasing signal; negative sensitivity represents movement away from it.
Biomass is clipped at zero after the update. The sensitivity has units
dependent on the signal concentration and spatial scale, so it must be
calibrated for each organism and field.

This is a spatial phenomenological capability, not a receptor-level model.
It does not yet represent receptor saturation, temporal adaptation, run-and-
tumble mechanics, or cell-level trajectories. Those should be added only when
the observation data can constrain them.
