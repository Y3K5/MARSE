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
| V2 | Resource-limited growth | Expected saturation or starvation behaviour | 3 | **Passing, well mixed**, in the batch engine and the version 2 network engine |
| V3 | Diffusion only | The numerical diffusion method, separately from microbial rules | 2 | **Passing, 1-D steady state** |
| V4 | Attachment and biofilm initiation | Transition from planktonic or seeded biomass to attached growth | 4 | Planned |
| V5 | Two-species competition | Expected dominance and coexistence regimes | 3 | **Passing, well mixed** |
| V6 | Cross-feeding | Explicit beneficial exchange shifts the equilibrium as designed | 3 | Planned |
| V7 | Environmental perturbation | Recovery or adaptation after a resource, pH or oxygen shift | 4 | Planned |
| V8 | Reproducibility | Re-running a saved manifest reproduces outputs within stated tolerances | 5 | **Passing** (batch, biofilm and ecosystem runs) |
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

## The two-dimensional ecosystem engine is not yet verified

Everything above concerns the well-mixed kernel and the one-dimensional biofilm
solver. The two-dimensional multispecies engine in `marse.ecosystem` is not
covered by any case in the table except V8: its runs now write a manifest and
`marse replay` reproduces them bit for bit, including a SHA-256 digest of the
whole final state. Reproducing a result is not the same as the result being
right, though. Its transport has not been checked against
`point_source_diffusion_2d`, and a review that ran the engine found the defects
below. Each
is written as a test in `tests/test_ecosystem_known_defects.py` that asserts
the *correct* behaviour and is marked as an expected failure. Fixing a defect
flips its test, which must then become an ordinary regression test.

| # | Defect | Measured | Consequence |
|---|---|---|---|
| 1 | The examples give oxygen a diffusivity of 10 µm²/h. The physical value in a biofilm is about 4×10⁶ µm²/h (`oxygen_diffusivity_um2_per_h` at 37 °C × 0.43). At 10 µm cells the explicit scheme could only run that value with 22 ms steps. | Oxygen spreads about 30 µm in 24 h instead of about 20 mm | Every oxygen gradient in the examples is set by the numerics, not the physics |
| 2 | Uptake follows *potential* growth, not actual growth. Production adds material that no consumed substrate pays for. | With growth held at zero, 22% of the carbon is consumed in an hour. With a yield of 1, 1.3–1.6 units are consumed per unit of biomass formed. | Yield and mass balance are both broken. Configuration schema version 2 already refuses such a process when it loads ([networks.md](networks.md)); the engine is fixed when it runs on version 2 The version 2 network engine meets the requirement: consumption is growth divided by yield, exactly, and nothing is consumed without growth |
| 3 | The periodontal species need oxygen to grow (a Monod term), although all three are anaerobes ([Holt and Ebersole 2005](https://pubmed.ncbi.nlm.nih.gov/15853938/)) | With no oxygen, growth is exactly zero | The study's oxygen-limited control points the wrong way |
| 4 | A species capability applies the substrate Monod term a second time | 50% of the intended rate at C = K, 13% at C = K/7 | Growth at low substrate is strongly understated. Version 2 refuses a second factor for one component, and a Monod factor gives exactly half the rate at C = K |
| 5 | Negative values are clipped instead of refused, and chemotaxis has no stability limit | One unstable chemotaxis step doubles total biomass | Mass is created silently. For reactions, the version 2 engine is positive without clipping at any step size; chemotaxis waits for the spatial engine |
| 6 | Mutation marks grid cells, including empty ones, not lineages, and every species draws from one shared random stream | At probability 1, every empty cell becomes "mutant" | Resistance does not move with the cells that carry it |
| 7 | Species are updated one after another within a step | Swapping two competitors in the file changes their final biomass by 0.9% | Results depend on how the file is written. In the version 2 engine every process acts on the same state: reordering changes fixed-step results only at rounding level |
| 8 | Each species has its own carrying capacity | Two species fill a cell to twice its capacity | Space is not shared. Shared space needs the spatial engine; a closed well-mixed box has none |
| 9 | ~~Ecosystem runs write no manifest and cannot be replayed~~ **Fixed:** they write a manifest (kind `ecosystem`) and replay exactly | — | The reproducibility claim now covers this engine; its manifests name the engine `ecosystem_v1_unverified` |

Until these are fixed, results from `marse ecosystem`, including the
periodontal study, support no quantitative or comparative conclusion. The
[roadmap](roadmap.md#order-of-work-correctness-first) fixes them before any
new mechanism is added.

The benchmark registry in `marse.validation.benchmarks` has one case scored
against simulator output so far: the first-order oxygen profile, met by the
validated one-dimensional solver in its first-order limit, at second order.
Its other cases have no simulator attached yet.

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

## Stoichiometric continuity

Reaction networks (configuration schema version 2) are checked when they load,
before anything runs: every process must conserve carbon, nitrogen and
electrons exactly ([theory.md §3.6](theory.md#36-composition-continuity-and-the-degree-of-reduction)).
The checks of that check, in `tests/test_schema_network.py` and
`tests/test_schema_formula.py`:

- **Textbook values.** Degrees of reduction of fifteen compounds, and the COD
  conversion factors 1.07 g g⁻¹ for glucose and 1.42 g g⁻¹ for cells.
- **Textbook reactions, exactly.** Respiration of glucose, alcoholic and
  homolactic fermentation, and nitrification, each derived from one given
  coefficient, water and protons included.
- **An independent method.** Growth rows equal the half-reaction method of
  Rittmann and McCarty (2001, ch. 2) in every coefficient, for three
  electron partitions. The transcribed half reactions are themselves checked
  for balance.
- **Atoms, not only the three quantities.** In forty random networks, every
  derived row balances carbon and nitrogen when recounted from the formulas,
  and oxygen, hydrogen and charge with the implicit water and protons.
- **What the engine will see.** Rounded once to double precision, every row
  still balances to within a few units in the last place. Reordering
  components, processes or `balanced_by` changes no coefficient.
- **Refusals.** Production from nothing (known defect 2's pattern), a row that
  nearly balances, and every malformed network are refused with a message
  naming the process, the field or the quantity.

## The version 2 network engine

A runnable network is integrated in a closed, well-mixed box
(`marse.core.well_mixed`; [theory.md §3.7 and §9.6](theory.md#37-rates)). The
checks, in `tests/test_well_mixed.py` and `tests/test_integrators.py`:

- **V2, against the analytical solution.** Monod batch growth reproduces
  `monod_batch_time` to a relative 1e-6. It ends at `batch_final_biomass`, and
  the integrator converges at second order, with the measured order above 1.9.
- **Balance over a long run.** Carbon, nitrogen and electrons stay within
  1e-12 of their totals over ten thousand steps. A planted leak of 1e-6 is
  caught at the first step, so the ledger is shown to be able to fail.
- **Positivity without clipping.** Steps a thousand times too long leave every
  concentration non-negative and every balance exact. The same holds for
  3,000 random networks with cycles, at steps up to 10⁶ h, with a worst drift of
  5e-16.
- **Selectivity.** A process that consumes nothing scarce runs at exactly the
  rate it would have alone.
- **The requirements behind known defects 2, 4, 5 and 7** hold on this engine,
  as the defect table above notes. The version 1 expected-failure tests stay
  until version 1 is removed.
- **Replay.** Runs replay bit for bit from their manifests. An edited manifest
  is refused.

## Running the suite

`python -m pytest` runs the pull-request suite, which leaves out tests marked
`slow` (long parameter sweeps). `python -m pytest -m slow` runs those, as the
nightly workflow does, and `python -m pytest -m ""` runs everything. As cases
land, `python -m pytest -m "numerical or invariance or regression"` runs the
fast checks and `python -m pytest -m benchmark` runs the biological
benchmarks. A `marse validate` command will run the same suite from an
installed package.
