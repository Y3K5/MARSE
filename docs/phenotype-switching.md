# Quorum-triggered phenotype switching

`PhenotypeConfig` defines a state activated by local biomass density relative
to carrying capacity by default. A species can instead declare an existing
`AdditiveConfig` as `quorum_signal` and a non-negative
`quorum_secretion_per_h`; biomass then secretes that signal, which is
transported and decayed as a continuous field. In signal mode, phenotype
thresholds use the signal field's declared concentration units.

Each state has an activation threshold, a lower or equal deactivation
threshold, optional growth and spreading multipliers, and a minimum dwell
time. The separate thresholds provide hysteresis so a noisy boundary does not
switch state every timestep.

The default implementation remains a phenomenological quorum-style mechanism.
Explicit secretion is a field-level source, not a receptor-level model. It
does not claim to model a specific signal molecule, regulatory network, or
intracellular protein abundance. Those details should be added only with
source-linked observations that can constrain them.
