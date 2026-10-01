# Uncertainty and sensitivity

`marse.analysis.uncertainty` provides a model-agnostic wrapper around deterministic
evaluators. A caller declares finite parameter ranges, supplies an evaluator
that returns one scalar outcome, and receives reproducible samples, output
quantiles, and rank-based screening correlations.

The default Latin-hypercube sampler gives each parameter coverage across its
declared interval while preserving a seeded random stream. `rank_sensitivity`
reports Spearman-style correlations, which are useful for monotonic screening
of growth, dose, transport, or ecosystem outputs.

These results are exploratory uncertainty propagation, not confidence or
credible intervals. They do not account for measurement noise, parameter
correlation, censoring, posterior distributions, or model-form uncertainty.
Those require an observation-specific statistical model and should not be
inferred from uniform ranges alone.
