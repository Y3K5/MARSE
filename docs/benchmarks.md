# Benchmark registry

`marse.validation.benchmarks` defines source-linked validation contracts.
Each case records its inputs, expected output, units, comparison metric,
tolerance, and source identifiers. A benchmark is a bounded claim about one
observable under one set of assumptions; passing it does not validate the
whole biological model.

The initial `analytical_benchmarks()` registry contains idealized exponential,
logistic, and first-order reaction-diffusion references. These cases protect
the numerical foundations while organism datasets are curated. Future
benchmarks can use the same comparison machinery for:

- strain-specific broth growth curves;
- CFU or colony-area trajectories;
- oxygen and pH profiles;
- additive or inhibition-zone responses;
- spatial ecosystem signatures.

Observed data, fitted parameters, and model outputs should remain separate.
Every experimental case must identify its source, endpoint units, replicate
structure, preprocessing, and tolerance rationale before it is promoted to
the registry.

The harness that runs a case and writes a scored report is described in
[phase2-plan.md](benchmarks/phase2-plan.md). It calls this registry when a
case supplies arrays. It does not replace it.
