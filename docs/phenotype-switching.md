# Quorum-triggered phenotype switching

`PhenotypeConfig` defines a state activated by local biomass density relative
to carrying capacity. Each state has an activation threshold, a lower or equal
deactivation threshold, optional growth and spreading multipliers, and a
minimum dwell time. The separate thresholds provide hysteresis so a noisy
boundary does not switch state every timestep.

The current implementation is a phenomenological quorum-style mechanism:
local biomass is the signal, and state changes are recorded per grid cell in
each frame. It does not claim to model a specific receptor, signal molecule,
regulatory network, or intracellular protein abundance. Those details should
be added only with source-linked observations that can constrain them.
