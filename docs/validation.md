# Validation plan

MARSE separates three questions. Is the software correct? Is the numerical
method correct? Does the model correspond to biology? Validation cases are
designed alongside the software, and benchmark definitions come before new
biology.

## Validation classes

| Class | Question | Example | Test marker |
|---|---|---|---|
| Unit correctness | Does each component implement its stated rule? | Growth update, neighbourhood lookup, boundary handling | *(none)* |
| Analytical or numerical | Does the method reproduce a known solution or convergence behaviour? | Diffusion against an analytical solution | `numerical` |
| Invariance | Are results unchanged under transformations that preserve the problem? | Grid translation or rotation where appropriate; seed control | `invariance` |
| Regression | Do changes unintentionally alter established results? | Golden benchmark trajectories | `regression` |
| Biological benchmark | Does the model reproduce selected observed behaviour? | Growth curve, response to nutrient limitation, biofilm maturation pattern | `benchmark` |
| Sensitivity | Which parameters drive a result? | One-at-a-time and global sensitivity analysis | |
| Uncertainty | How stable are conclusions across plausible parameter ranges? | Ensembles with confidence intervals | |

## Anatomy of a validation case

Every `ValidationCase` records:

- **ID and question:** the behaviour being checked.
- **Setup:** a small, synthetic experiment configuration committed with the case.
- **Target metric:** what is measured, with units.
- **Tolerance:** the acceptable deviation, and why.
- **Reference:** the analytical solution, dataset or publication it is compared against.
- **Status:** planned, implemented or passing.

## First benchmark suite

| ID | Case | What it demonstrates | Phase | Status |
|---|---|---|---|---|
| V1 | Single-species unrestricted growth | The chosen growth law and carrying-capacity behaviour | 3 | **Passing, well mixed** |
| V2 | Resource-limited growth | Expected saturation or starvation behaviour | 3 | **Passing, well mixed** |
| V3 | Diffusion only | The numerical diffusion method, separately from microbial rules | 2 | **Passing, 1-D steady state** |
| V4 | Attachment and biofilm initiation | Transition from planktonic or seeded biomass to attached growth | 4 | Planned |
| V5 | Two-species competition | Expected dominance and coexistence regimes | 3 | **Passing, well mixed** |
| V6 | Cross-feeding | Explicit beneficial exchange shifts the equilibrium as designed | 3 | Planned |
| V7 | Environmental perturbation | Recovery or adaptation after a resource, pH or oxygen shift | 4 | Planned |
| V8 | Reproducibility | Re-running a saved manifest reproduces outputs within stated tolerances | 5 | **Passing** |
| V9 | IWA benchmark BM1 | Substrate flux and concentration for a monospecies biofilm at fixed biomass, against published reference solutions | 3 | **Now reachable** |
| V10 | IWA benchmark BM3 | Multispecies, multisubstrate biofilm (heterotrophs, nitrifiers, inert biomass) | 4 | Proposed |

**"Well mixed" is a real qualifier, not a hedge.** V1, V2 and V5 are verified
for a homogeneous batch or chemostat, where the analytical solutions apply.
Their spatial counterparts — the same behaviours inside a biofilm with
gradients — are what Phases 2 and 3 add, and they can diverge from the
well-mixed answer substantially. V5 is the clearest example: competitive
exclusion is a well-mixed result, and spatial structure can permit the
coexistence it forbids (theory.md §7.3).

V3 covers the one-dimensional steady-state solver
(`marse.spatial.diffusion`). It is checked against three analytical limits
that bracket real Monod uptake — no uptake, first order and zero order — plus
a flux balance and a convergence study confirming the discretisation is
second order. The transient and two-dimensional cases, for which
`point_source_diffusion_2d` is the reference, arrive with biofilm structure.

One note from building it, because it generalises. Comparing a numerical
solution against an *approximate* analytical solution carries two errors: the
discretisation error, which shrinks under refinement, and the error in the
approximation itself, which does not. In the first-order limit the second is
of order `C / K`, and at `K = 1e4` it sat above the grid error and flattened
the measured convergence to zero — the solver looked first-order-at-best when
it is second order. Any convergence study here must first establish that the
reference is closer to exact than the grid error being measured.

V9 and V10 adopt the published benchmark problems of the IWA Task Group on
Biofilm Modeling, which exist precisely to compare modelling approaches and
have published reference solutions from multiple independent models. Using
them is stronger evidence than any benchmark MARSE could define for itself.
See [modeling-landscape.md §5.1](modeling-landscape.md#51-the-iwa-benchmark-problems).

## Structural output metrics

Simulated biofilm structure is reported using the vocabulary established for
confocal microscopy analysis, so that simulations and images are directly
comparable: **biovolume, mean and maximum thickness, substratum coverage,
roughness coefficient** and **volume-to-surface ratio**
([modeling-landscape.md §5.2](modeling-landscape.md#52-structural-metrics)).
These give validation case V4 a quantitative target it currently lacks.

## Stratified growth

Coupling the solute field to the biomass that consumes it
(`marse.biofilm.biomass`) reproduces the behaviour that motivates spatial
modelling in the first place. With identical cells at uniform density, growth
is confined to a surface layer whose depth stops increasing once the biofilm
passes the penetration depth, so the biofilm-averaged growth rate falls
roughly as 1/thickness while the surface keeps growing at its unlimited rate:

| Thickness | Active zone | Mean growth rate | vs. surface |
|---|---|---|---|
| 25 µm | 25 µm | 0.299 h⁻¹ | 100% |
| 100 µm | 100 µm | 0.298 h⁻¹ | 99.7% |
| 200 µm | 137 µm | 0.200 h⁻¹ | 66.9% |
| 400 µm | 137 µm | 0.100 h⁻¹ | 33.4% |
| 800 µm | 137 µm | 0.050 h⁻¹ | 16.7% |

Two limits are checked rather than just the interesting one: a film thinner
than the penetration depth must reduce to the well-mixed answer, and a thicker
one must stratify. This is qualitatively the narrow band of protein synthesis
reported for real biofilms; turning it into a quantitative benchmark is what
V9 (IWA BM1) is for, and BM1's fixed-biomass monospecies setup is exactly the
problem this code now solves.

## Analytical reference solutions

Cases V1, V2, V3 and V5 have exact solutions, implemented in
`marse.validation.analytical` and derived in [theory.md](theory.md):

| Reference | Problem | Theory |
|---|---|---|
| `exponential_growth` | Unrestricted growth | §1.1 |
| `logistic_growth` | Growth to a carrying capacity | §1.2 |
| `monod_batch_time` | Integrated Monod batch culture | §3.5 |
| `batch_final_biomass` | Final biomass from mass balance | §3.4 |
| `zero_order_penetration_depth` | Solute penetration into a biofilm | §5.1 |
| `point_source_diffusion_2d` | Diffusion from a point release | §4.2 |
| `chemostat_break_even` | Break-even substrate level, the $R^{*}$ rule | §7.1 |

A wrong reference would silently validate a wrong simulation, so each one is
itself checked against an independent numerical solution — fourth-order
Runge–Kutta for the growth problems, a Newton-solved nonlinear
reaction–diffusion system for the penetration depth, finite differences for
diffusion, and a full three-species ODE integration for the chemostat. These
checks run as part of the ordinary test suite.

## Running the suite

As cases land, `python -m pytest -m "numerical or invariance or regression"`
runs the fast checks and `python -m pytest -m benchmark` runs the biological
benchmarks. A `marse validate` command will run the same suite from an
installed package.
