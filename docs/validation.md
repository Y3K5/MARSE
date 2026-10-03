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
| 1 | The examples give oxygen a diffusivity of 10 µm²/h. The physical value in a biofilm is about 4×10⁶ µm²/h (`oxygen_diffusivity_um2_per_h` at 37 °C × 0.43). At 10 µm cells the explicit scheme could only run that value with 22 ms steps. | Oxygen spreads about 30 µm in 24 h instead of about 20 mm | Every oxygen gradient in the examples is set by the numerics, not the physics. The version 2 spatial engine runs oxygen at its physical diffusivity: transport is implicit, so no step limit applies ([theory.md §9.8](theory.md#98-implicit-reactiontransport-integration)) |
| 2 | Uptake follows *potential* growth, not actual growth. Production adds material that no consumed substrate pays for. | With growth held at zero, 22% of the carbon is consumed in an hour. With a yield of 1, 1.3–1.6 units are consumed per unit of biomass formed. | Yield and mass balance are both broken. Configuration schema version 2 refuses such a process when it loads ([networks.md](networks.md)), and the version 2 network engine meets the requirement: consumption is growth divided by yield, exactly, and nothing is consumed without growth |
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

## Transport in one, two and three dimensions

The spatial engine moves dissolved components between cubic voxels by
finite-volume diffusion (`marse.spatial.transport`,
[theory.md §4.7](theory.md#47-finite-volume-transport-on-voxels)). The checks,
in `tests/test_transport.py`:

- **Exact eigenvectors.** With lateral faces periodic, no flux through the
  substratum and the bulk held at zero, a product of cosines is an exact
  eigenvector of the discrete operator in 1-, 2- and 3-D. Its rate matches
  `cosine_mode_rate` to rounding, and that rate converges to the continuum
  rate at second order.
- **The 3-D point release.** On a spreading Gaussian
  (`point_source_diffusion_3d`), the operator's error against the exact time
  derivative falls at second order: the measured orders are 1.93 and 1.98 as
  the voxels shrink from 4 µm to 1 µm.
- **Conservation.** The rate summed over the box equals the flux through the
  top face to rounding, for random fields in every dimension. A component
  without diffusivity does not move.
- **Symmetry and dimension.** A laterally uniform 3-D field diffuses exactly,
  bit for bit, as a single column. Swapping or reflecting the lateral axes
  moves the result with them.
- **The multigrid solver** (`marse.spatial.multigrid`).
  - Each V-cycle reduces the error by less than 0.15, the same on grids from
    8³ to 32×32×16 and in 1-D and 2-D. On 64×64×32 the measured factor is
    0.06–0.09.
  - GMRES with this preconditioner agrees with a direct solve.
  - The same system gives the same answer, bit for bit.
  - The coarsest level's inverse, from MARSE's own blocked LU factorisation,
    agrees with LAPACK's to 1e-12, including a matrix that needs pivoting. A
    singular matrix is refused. Its dense matrix, assembled with one
    application per voxel, equals the operator applied to each unknown in
    turn, bit for bit.
- **Output.** VTK frames read back exactly in double precision and to 1e-7 in
  single precision.

## Reactions and transport in space

The version 2 spatial engine (`marse.core.reactive_transport`,
[theory.md §9.8](theory.md#98-implicit-reactiontransport-integration)) runs a
network in a box of voxels over a surface. It was prototyped first in one
dimension against criteria set in advance, and it passed all of them:

| Criterion | Measured |
|---|---|
| C, N and e⁻ conserved over 10⁴ steps, counting imports | to 2.4e-15 of the totals |
| Nothing negative; the limiter rarely needed | never negative; the limiter was used in 1 of 10⁴ steps |
| Second order in time | order 1.74 → 1.92 as the step falls from 90 s to 2.8 s (stiff order reduction) |
| A sudden drop in bulk oxygen followed | 90 s steps within 3e-6 of air saturation |
| Long steps reach the quasi-steady profile | within 7e-7, after six minutes, of the same grid's steady state |

The tests, in `tests/test_reactive_transport.py`:

- **The scheme's own algebra.** On a diffusion mode, each step multiplies the
  mode by exactly the stability function of the method, to 1e-9, and the
  result converges to the exact decay at second order.
- **Agreement with the well-mixed engine.** Without diffusion, every voxel of
  the box evolves as the well-mixed engine (2b) evolves the same contents, to
  1e-5.
- **Agreement across dimensions.** A laterally uniform 3-D box, taking the
  same steps as a single column, evolves as that column to 1e-10 of each
  component's peak (measured: 2e-13, rounding), and no lateral flux arises.
- **The ledger in an open box.** Carbon, nitrogen and electrons equal their
  imports through the top to 1e-12 of the totals. A planted leak in the
  transport is caught at the first step, so the ledger is shown to be able to
  fail.
- **Positivity.** Hour-long steps over a biofilm with an anoxic base leave
  every concentration non-negative. Nothing is clipped: a ledger held to
  1e-12 balances after every step. A trace of four units in the last place,
  sending one through each of six faces, cannot be scaled into balance,
  because each scaled transfer rounds back up to one unit. It stays
  non-negative with its total unchanged, both in the limiter's rounds and in
  its fallback.
- **The analytic Jacobian** of the rates matches finite differences of the
  rates as the engine evaluates them, including below zero, where it has no
  slope.
- **The same result on any number of threads.** The 1-D example, run in
  subprocesses with OpenBLAS at 1 and at 4 threads, ends with one final
  digest. Before the coarsest level left LAPACK it did not: the digests
  differed at 1, 2 and 4 threads, so a run replayed only on a machine with as
  many cores.
- **Schema, replay and command line.**
  - Every impossible domain is refused with a message that says why.
  - Random colonies are placed by the seed and replaced by it.
  - Runs replay bit for bit, and an edited manifest is refused.
  - `marse run` writes the totals, the frames and the ParaView files.

### Performance

The stage pre-registered a target: a 64 × 64 × 32 box (131,072 voxels,
128 × 128 × 64 µm) should run 24 simulated hours in under an hour on one core.
**It missed: the run took 2 h 17 min (8,227 s).**

- **The scene:** the 3-D example with 2 µm voxels (`voxels` [64, 64, 32],
  `voxel_um` 2, `duration_h` 24).
- **The machine:** one core of a 2.8 GHz Intel Xeon cloud container, with
  NumPy's BLAS held to one thread.

Where the time went:

- **The start-up, 41 minutes (about 2,440 s) for the first 15 simulated
  minutes.** The colonies start in fresh liquid, and the start-up transient
  takes several hundred steps at about 7.7 s each. Up to 51 of the limiter's
  rounds per step went to traces far below any tolerance.
- **The remaining 23.75 hours, 5,786 s.** That is 244 s per simulated hour on
  average, falling from 330 s early on to 156 s over the last eight hours as
  the steps lengthened.
- **In all, 691 steps**: 3 rejected, 105 limited, no Newton failures, and 4,324
  Newton iterations.
- **Conservation held throughout:** carbon, nitrogen and electrons balanced to
  3.8e-16.
- **Scaling:** a step costs 16 times as much as on the 32 × 32 × 16 example,
  for 8 times the voxels.

By 24 h the biomass has grown in place to about 55 times its starting amount,
which is not physical until colonies spread (Stage 2d). The run measures the
solver, not the biology. Stage 7 takes up the acceleration of the start-up and
of the cost per voxel, without loosening any tolerance
([roadmap](roadmap.md)).

## Adhesion to surfaces

Cells in the liquid bind to the substratum, detach from it or lock onto it
([theory.md §6.4](theory.md#64-attachment-to-surfaces)). The checks, in
`tests/test_surfaces.py` and `tests/test_adhesion.py`, were set before the
code was written:

| Criterion | Threshold | Measured |
|---|---|---|
| The Lévêque constant, and the wall flux of the boundary-layer problem marched numerically | 1e-12; 1% | exact to rounding; within 1% |
| Stokes–Einstein, a 1 µm sphere in water at 25 °C | 0.4907 µm² s⁻¹ to 0.1% | 0.4907 |
| The engine against the closed-form kinetics with blocking, in 1-, 2- and 3-D | the solver's tolerance | 2.0e-6 relative at a tolerance of 1e-6, identical in every dimension |
| The Langmuir limit without locking, $n_\infty = j_0/(j_0/n_J + k_{\mathrm{off}})$ | 1e-6 | 7e-16 |
| Coverage never passes the jamming limit, and reaches it without detachment | exact; 1e-3 | never above 0.547; reaches it to 1e-3 |
| Each patch of a patterned box evolves as its own column | rounding, 1e-12 | agrees to 1e-12 |
| The ledger, with deposition and detachment | 1e-12 of the totals | 4e-16 to 6e-16 |
| Hour-long steps | never negative, nothing clipped | never negative; the ledger balances after every step |
| Without a substratum, a run of Stage 2c | bit for bit | bit for bit: the 1-D and 3-D examples end with the same digests with and without E1's code |
| Replay, and the schema's refusals | bit for bit, on any number of threads; each refusal says why | as required; both scenes end with one digest at 1 and 4 threads |
| The lab scene's initial deposition | inside 0–2.9 × 10³ cm⁻² s⁻¹ (Sjollema et al. 1988) | 595 on bare glass, 357 on coated glass |
| The dental scene, zirconia over titanium | 0.63 ± 1% after 2 h; 0.60–0.72 covered after 24 h | 0.63; 0.67 |

The last row checks consistency with a calibration, not a prediction:
zirconia's attachment efficiency was set to 0.63 of titanium's from Scarano et
al. (2004) ([parameters.md §7.2](parameters.md#72-titanium-against-zirconia)).
The covered ratio drifts to 0.67 over the day because the covered area is not
proportional to the cells bound. Locked cells grow, and a surface grown past
jamming binds no more.

The tests in detail:

- **Delivery** (`marse.spatial.colloids`). A cell's delivery is checked three
  ways:
  - the Lévêque constant is checked against its closed form;
  - the wall flux of the convection–diffusion equation, marched downstream on
    a stretched grid, matches it to 1%;
  - the transfer velocity scales as $D^{2/3}$, $\dot\gamma^{1/3}$ and
    $x^{-1/3}$.

  Stokes–Einstein matches the known diffusivity of a micron sphere, and the
  salivary film shears its surface at 0.4 to 5.4 s⁻¹.
- **The substratum** (`marse.spatial.surface`). Patches give every face
  exactly one material, by where its centre lies. A pattern with a gap, an
  overlap, a patch outside the box or covering no face, or the wrong number of
  bounds is refused with the reason.
- **The closed form.** `adhesion_kinetics` matches an independent fourth-order
  Runge–Kutta integration, with and without locking and from a surface
  already partly covered. It approaches the jamming limit and never passes it.
- **The engine against the closed form.** A box with no dissolved components
  and no growth binds cells as `adhesion_kinetics` says, in 1-, 2- and 3-D.
  Without locking it settles on the Langmuir balance, and with it every
  surface jams.
- **Geometry.** A box patterned with two materials, taking the same fixed
  steps as two single columns, evolves face by face as those columns do.
- **Conservation.** Carbon, nitrogen and electrons equal what crossed the top
  and the substratum, to 1e-12 of the totals. With nothing growing, the cells
  on the surface equal the cells imported through it.
- **Positivity.** Hour-long steps leave every concentration non-negative, and
  a ledger held to 1e-12 balances after each one.
- **The analytic Jacobian** of the exchange and of locking matches central
  finite differences to 1e-6. The dental scene is sampled at random, with 20%
  of the values negative and faces on both sides of jamming.
- **The schema.** A domain without a surface writes exactly the fields it wrote
  before, so its run id and checksums are unchanged. One with a surface
  round-trips. Thirteen impossible scenes are refused with the reason:
  - a missing field;
  - one component named as both states of a species, or two states with
    different formulas;
  - a reversible state that diffuses;
  - an unknown material, or a species and material pair stated twice or not
    at all;
  - an efficiency outside [0, 1], or a negative rate;
  - a temperature below absolute zero, or a shear or cell diameter of zero;
  - patches that overlap.

  A reversible state that takes part in a process is refused too.
- **The scenes and the command line.** The lab scene binds at an initial rate
  inside the published range. The dental scene binds to zirconia at 0.63 of
  titanium, and to enamel exactly as to titanium. `marse run` writes
  `surface.csv` and replays it bit for bit. `marse check` previews delivery
  and binding on every material. A day on the dental surfaces, which covers
  titanium more than zirconia, runs in the slow suite.

### Timings

On one machine (4 cores, Python 3.12):

| Run | Voxels | Simulated | Wall time |
|---|---|---|---|
| Lab flow chamber, `marse run` | 8 × 4 | 4 h | 0.5 s |
| Dental surfaces, `marse run` | 16 × 4 × 8 | 24 h | 35 s |
| Dental surfaces as four columns, `examples/surface_adhesion.py` | 8 per column | 24 h each | 6.6 s |

The times are the same at 1 and 4 BLAS threads, to within a few percent, and
so are the results, bit for bit.

## pH from electroneutrality

A network whose components include acids and bases sets a pH in every voxel,
from the charge balance over them
([theory.md §3.8](theory.md#38-acidbase-equilibria-and-ph)). These criteria
were set before Stage S1 was built, and are checked in `tests/test_acid_base.py`:

| | Criterion | Threshold | Result |
|---|---|---|---|
| G1 | The hydrogen ion concentration against an independent bisection of the charge balance, for lactate, carbonate, phosphate and ammonium, 500 random compositions each | 1e-10 relative | 2.4e-12 |
| G2 | The charge left in 2,000 random plaque-like voxels, pH 1 to 13 | 1e-12 of the potassium, chloride and hydrogen ions | 8.1e-14 |
| G4 | The rate Jacobian of a process with a pH factor, through dh/dc, against finite differences | 1e-6 relative | 2.1e-7 |

Other checks:

- **Closed forms.** A weak acid, with the potassium that neutralises it at a
  chosen pH, gives back that pH to 1e-12, from pH 3 to 8.5. A half-neutralised
  acid sits at its pKa. Pure water sits at half of pKw.
- **The slopes.** dh/dc agrees with Richardson-extrapolated finite
  differences to 1e-7. Acids and anions raise h, cations lower it, and
  uncharged components do not move it.
- **Convergence to the last bit.** The solve stops where the balance reaches
  its own rounding floor. Without that test, Newton's method stepped for ever
  between two neighbouring values of h in a lactate buffer at pH 7.
- **No cycles.** In one plaque voxel of an oral scene, every Newton step
  landed on the other end of an unchanging bracket, between pH 3.00 and 5.76,
  and the solve never finished. A step must now halve the one before, or the
  bracket is bisected (Press et al. 2007). That voxel converges to the root
  bisection finds, and 2,000 random plaque voxels, pH 1.3 to 12.7, take 12
  iterations, against 19 without the safeguard.
- **The dissociated factor**, the protons a total has lost at the local pH,
  agrees with its closed forms for one pKa and for phosphate, and its slope
  with finite differences in log h. The rate Jacobian through it agrees with
  finite differences to 5.9e-8. At equilibrium, the cations the fixed groups
  hold are their charge.
- **The cardinal pH slope** agrees with finite differences of the factor.
- **Runs.**
  - A fermenting well-mixed box acidifies at every record and writes its pH.
  - A plaque column writes `ph.csv`, from pH 7.00 at the start.
  - `marse check` prints the pH each starting composition implies.
- **The elements.** Phosphate and the salt ions hold no electrons. A process
  that makes potassium from chloride is refused. `balanced_by` closes an
  element.
- **Nothing else changes.** The examples' final digests are bit-identical to
  those before pH existed, and one absolute tolerance given per component
  reproduces the run with a single tolerance, bit for bit, in a box and in
  space.

## The mouth and its film

A domain under a salivary film exchanges with the mouth's saliva, whose
composition is solved with the box
([theory.md §4.8](theory.md#48-a-salivary-film-and-the-mouth) and
[§9.9](theory.md#99-the-mouth-solved-with-the-box)). The criteria set before
Stage S1 was built, checked on a column of 300 µm of plaque under 100 µm of
film, in `tests/test_mouth.py`:

| | Criterion | Threshold | Result |
|---|---|---|---|
| G6a | With a constant flow and nothing taken up, the mouth's sugar after each swallow against Dawes's closed form, (H_resid / H_max) to the number of swallows | rounding | 3.7e-11 over ten swallows |
| G3, G6b | The box and the mouth together conserve every quantity over an hour with a residue of 10% sucrose, counting what was secreted and swallowed | 1e-12 | 5.6e-16; the box alone 3.0e-16 |
| G6c | The pH at the substratum over that hour, with the mouth running 60 s and 5 s ahead of the box | 1e-3 | 6.1e-5 |

That hour takes 5.9 s: 358 steps and 71 swallows.

Other checks, in `tests/test_reservoir.py`, `tests/test_mouth.py` and
`tests/test_transport.py`:

- **A film and its pool** relax towards their volume-weighted mean at
  k (1 + delta / H), as the closed form says, to 1e-6.
- **The bordered linear system** is the Jacobian of the step, by central
  differences, in one and two dimensions. It is solved exactly.
- **The pool never goes below zero.** Asked for more than it holds, it gives
  what it has, and both ledgers still balance.
- **At rest**, a column of plaque under resting saliva stays at pH 7.00 and
  swallows once a minute.
- **The film** is renewed fastest at its surface, not at all below it, and at
  u_bar / l on average.
- **A closed top** keeps everything in the box. A column's direct solve agrees
  with multigrid, and the operator it applies is multigrid's, bit for bit.
- **Replay** reproduces a run under the mouth exactly, and every impossible
  film or mouth is refused with the reason.

## The diet

A mouth takes the rinses, drinks and foods of a diet, and the food some of
them leave on the teeth ([theory.md §4.9](theory.md#49-the-diet)). On the
column of the section above, in `tests/test_mouth.py`:

| | Criterion | Threshold | Result |
|---|---|---|---|
| G3, G6b | The box and the mouth together conserve every quantity over an hour with a rinse of 10% sucrose, 200 mL of a sugared drink sipped over 20 minutes and a sweet that leaves food on the teeth, counting what was secreted, eaten, swallowed and expelled | 1e-12 | 1.2e-15; the box alone 8.6e-16 |
| G6c | The pH at the substratum over an hour after a rinse of 10% sucrose held for a minute, which mixes the film with the mouth once a second, with the mouth running 60 s and 5 s ahead of the box | 1e-3 | 8.6e-5 |

Over the same hour, the mouth's sugar agrees within 8.4e-4. Without the
run-ahead's expectation of what the plaque gives back (theory.md §9.9), it
agreed within 1.2e-2 and the pH within 2.3e-4. The hour after the rinse takes
6.5 s: 404 steps and 73 swallows. The hour with the drink takes about 30 s,
because a drink sipped at 10 mL a minute makes the mouth swallow 910 times,
and every swallow ends a span.

Other checks:

- **What each intake brings** is booked as eaten exactly: the sugar of the
  rinse, the drink and the sweet, the drink's salts, and the food left on the
  teeth, within 5e-15 of what was stated.
- **A rinse** is held without a swallow while the glands secrete into it,
  and is expelled to the resting volume at its end, keeping the resting
  share of every amount.
- **A drink** makes the mouth swallow for every 0.3 mL it and the saliva
  add, and holds the mouth's sugar at q c / (q + Q) of the drink's, within 2%.
- **A sweet** releases what it states, adds no liquid, and raises the flow
  above 1 mL a minute as it is tasted.
- **Food left on the teeth** goes into the film over its region only, in two
  dimensions, and dissolves as the network's process releases it.
- **Mixing** during a rinse brings more of its sugar into the plaque than the
  film's renewal alone.
- **A diet that starts after the run ends** changes nothing: the final state
  is bit for bit that of the run without it.
- **Replay** reproduces a run with a diet exactly, `marse check` lists the
  diet, and every impossible diet is refused with the reason.

## The Stephan curve

The oral scenes of `examples/environments/oral`
([environments.md](environments.md#the-oral-scenes)): 150 µm of plaque under
a salivary film and a mouth, given a rinse of 10% sucrose, 100 mL of a drink
of 10% sucrose sipped over 20 minutes, or the rinse and food left on the
teeth. The last two criteria set before Stage S1 was built, in
`tests/test_oral_scenes.py`, as slow tests:

| | Criterion | Threshold | Result |
|---|---|---|---|
| G5 | After the rinse, the pH at the substratum falls at least one unit, to a minimum within 5 to 20 minutes, is back above 6 within an hour, and the plaque holds more lactate at 7 minutes | 1 unit; 4.5 to 5.5; 60 min; 10 to 60 mM | 1.93, from 6.80 to 4.87 at 15.6 min; back above 6 at 42.9 min; +15.0 mM |
| G7 | A two-hour curve in 100 voxels | under 30 s | 12.7 s: 550 steps and 132 swallows |

Both ledgers hold to 3.1e-15 or better in every scene, and the film is never
more acid than the plaque below it.

**What calibration took.** The prototype met G5 with a plaque buffer whose
cations were held fixed. In MARSE, the scenes first did too, but the film
above the plaque then fell to pH 3.3: the lactate leaving the plaque carried
the acid's protons out with it. With the buffer's cations released as its
groups take up protons (theory.md §3.8), the acid stays in the plaque until
the saliva's alkali takes it off, and 300 µm of plaque stayed acid for more
than an hour. Four values were calibrated again, within their ranges: the
plaque is 150 µm thick, the film moves at 6 mm per minute, the buffer is
160 mM, and every charged component diffuses at 7 × 10⁻¹⁰ m² per s
([parameters.md §8](parameters.md#8-saliva-plaque-and-diet)). All four are
confidence C.

**Directions**, each against the same rinse over 90 minutes, each a slow
test:

| Change | Lowest pH | Minutes below pH 5.5 | Back above pH 6 | As reported |
|---|---|---|---|---|
| None: the rinse | 4.87 at 15.6 min | 31.7 | 42.9 min | the Stephan curve (Stephan 1944) |
| 100 mL sipped over 20 min | 4.80 at 34.8 min | 52.3 | 64.7 min | a continuous presence of sugar keeps plaque acid |
| Food left on the teeth, 0.02 mol/m² | 4.78 at 22.2 min | 56.2 | 69.5 min | retained food keeps it acid (Kashket, Zhang and Van Houte 1996) |
| Low flow: 0.1 mL/min, up to 1.0 more | 4.84 at 18.0 min | 39.7 | 52.0 min | a greater fall, low for longer (Lingström and Birkhed 1993) |
| A slower film: 2 mm/min | 4.42 at 33.0 min | 87.1 | not within 90 min | a lower minimum, a slower return (Macpherson and Dawes 1991) |
| Thicker plaque: 300 µm | 4.57 at 22.2 min | 85.7 | not within 90 min | acid for longer (Dawes and Dibdin 1986) |
| More fixed buffer: 240 mM | 5.03 at 15.6 min | 31.0 | 47.3 min | a shallower fall, a longer low phase (Dibdin 1990) |

Dawes and Dibdin (1986) found an optimum thickness at which plaque reaches
its lowest pH, so a thicker plaque need not fall further; the test asks only
that it stays acid longer. `python examples/stephan_curve.py` runs the first
three and checks G5 and the first two directions (44 s).

## Plaque that spreads

Plaque in a column spreads up it as it grows, is worn at its surface and
detached above a maximum height, and is brushed off
([theory.md §6.1](theory.md#61-biomass-balance) and
[§6.3](theory.md#63-detachment)). The criteria set before Stage S2 was built,
in `tests/test_spreading.py`:

| | Criterion | Threshold | Result |
|---|---|---|---|
| P2 | After every packing, no voxel is filled past full, and the front is sharp: full voxels, then at most one partly filled. Checked on 500 random columns, against the remap done voxel by voxel | rounding | at most 4.7e-15 past full; at most one partly filled voxel |
| P3 | A film growing at 0.1 per hour and worn at 4 µm per hour, from 20 µm and from 60 µm, follows L(t) = u/μ + (L₀ − u/μ) e^{μt} (Wanner and Gujer 1986) | 1% | 9.7e-7 and 3.1e-8 |
| P4 | Under the mouth, a brushing 15 minutes into a half hour takes 42% of 150 µm of plaque off from the surface down, and the mouth expels it | rounding | 87.000000000 µm left; what was removed within 1.1e-16 of 42% |
| P1, in part | The box and the mouth together conserve every quantity over that half hour, counting what was brushed off and expelled | 1e-12 | 1.9e-16; the box alone 5.1e-16 |

The half hour takes 1.1 s: 30 swallows. P1 in full, with oxygen and chewing,
comes with them.

Other checks:

- **Wear alone** takes the surface off at its velocity, to 1.4e-14 µm, and
  books all of it as detached, with the buffer's groups it carried.
- **Growth past the maximum height** is detached. Over six hours, 1.0e-3 more
  was detached than in the continuous solution, the first-order error of
  cutting after each step (theory.md §6.1).
- **The splitting of wear** matters. With half the wear before each step and
  half after, P3's error follows the integrator's tolerance. With all of it
  after, the error is 3.6e-3 at the same tolerance (theory.md §6.3).
- **Carried components**, such as the buffer's groups, move with the solid
  and are never lost. Where a voxel holds no solid, they stay.
- **Food left on the teeth** loses the same share to a brushing as the
  plaque, to 1e-12, and the mouth expels it. A plaque that lists that food
  as part of itself is refused.
- **A plaque that grows nothing** stays packed where it was, so S1's rinse
  runs as before, within 1e-12.
- **The books** close to rounding in every case: every ledger to 1.2e-15 or
  better.
- **Replay** reproduces a run whose plaque spreads exactly. `marse check`
  describes the plaque and its cleanings, and every impossible plaque or
  cleaning is refused with the reason.

## Oxygen from the air

The top face at the air, which holds each gas it lists at its saturation
([theory.md §4.10](theory.md#410-oxygen-from-the-air)). The criterion set
before Stage S2 was built, in `tests/test_oxygen.py`:

| | Criterion | Threshold | Result |
|---|---|---|---|
| P5a | Oxygen filling a slab 200 µm deep from the air, against the series solution (Crank 1975), at 2, 7, 18 and 72 s | 2% after the grid is refined | at most 2.2e-4 of saturation on 5 µm voxels and 5.5e-5 on 2.5 µm: second order, a ratio of 4.06 |
| P5b | Under uptake of zero order, oxygen reaches δ = √(2DC_s/k₀) = 150 µm, falling to 1% of saturation 135.0 µm in | 2% | 135.3 µm on 5 µm voxels; the profile within 2.7e-4 of saturation of (1 − x/δ)², the difference being the Monod saturation constant of 1e-5 mol m⁻³ near the front |

Both ledgers count what the air gave and took. In the column, the air's
book and the oxygen gained agree to 4.7e-17. Under the mouth, S1's rinse with
plaque that respires sugar closes the box to 9.5e-16 and the box and the
mouth together to 3.5e-16, over 15 minutes and 27 swallows; the mouth's
oxygen stays at saturation throughout, whatever the swallows and the rinse
take or bring. The film's surface sits within 0.3% of saturation, and the
oxygen falls from there into the plaque.

Other checks:

- **Nothing else crosses the face.** A component the air does not list
  neither enters nor leaves through it.
- **A rinse that brings oxygen of its own** is brought back to saturation,
  and the difference is booked.
- **What is refused:** a gas that is particulate, an ion or an acid-base
  total (carbon dioxide stays with carbonate), a negative saturation, a gas
  also stated in the saliva, and a bulk liquid at the air.
- **Replay** reproduces a run at the air exactly, and `marse check` describes
  the air.

## Chewing

A chewed food adds the mouth's chewing flow while it lasts
([theory.md §4.9](theory.md#49-the-diet)). In `tests/test_mouth.py`:

- **The flow** is 1.3 mL per minute while 1.0 mL per minute of chewing is
  added to the resting 0.3, and back to 0.3 when it ends, to rounding. The
  mouth swallows 43 times in the 10 minutes of chewing, and 5 times in the 5
  minutes after.
- **The saliva** comes half way to stimulated saliva at that flow: the
  mouth's carbonate rises from 5.0 to 10.0 mM.
- **The books** close to 5.0e-16 for the box and 1.8e-16 for the box and the
  mouth together.

The criterion set before Stage S2 was built, as a slow test on S1's rinse
with 1 mL per minute of chewing:

| | Criterion | Threshold | Result |
|---|---|---|---|
| P7 | Sugar-free gum chewed after a rinse of 10% sucrose brings the plaque's pH back sooner (Dibdin, Dawes and Macpherson 1995) | the direction | Chewed from minute 2, the fall stops at pH 5.64 at 2.4 minutes, and the plaque is back above 6 at 3.9 minutes. Chewed from minute 10, it is back above 6 at 13.7 minutes. The rinse alone: 4.87 at 15.6 minutes, back at 42.9. |

The size of the effect rests on mixing the film once a second while gum is
chewed, as for any intake, which is confidence C.

## A day of plaque

The two scenes Stage S2 added to `examples/environments/oral`
([environments.md](environments.md#a-day-of-plaque)): oxygen in 400 µm of
plaque before and after a rinse of 10% sucrose, and a day of meals, sweets,
gum, two brushings, wear and the air. The last criteria set before Stage S2
was built, in `tests/test_plaque_scenes.py`, as slow tests:

| | Criterion | Threshold | Result |
|---|---|---|---|
| P6 | Under saliva, plaque is anoxic below about 220 µm; after sucrose, oxygen reaches less far, about 150 µm (von Ohle et al. 2010) | 200 to 250 µm; shallower after sucrose | 216 µm under saliva; 151 µm, 9 minutes after the sucrose; back to 216 µm 45 minutes later |
| P1 | The box, and the box and the mouth together, conserve every element over a day with growth, packing, wear, brushing, chewing and the air | 1e-12 | 9.8e-15 for the box; 4.6e-16 for the box and the mouth |
| P8 | The cost of that day in a column of 100 voxels | set from S1's rate, 152 s a day | 130 s: 3,000 steps, 2,156 swallows |

**What calibration took.** In the prototype, the plaque respired its own
biomass between meals. Over a day that thinned it from 150 to 33 µm, since it
grew only on meals. In MARSE it lives on saliva instead. Saliva's mucins carry
about 1.5 mM of hexose (Payment et al. 2000; Levine et al. 1987), of which the
scenes take 1 mM as available. Where oxygen reaches, the plaque grows on it
with a yield of 3 C-mol per hexose, and respires a little of it for
maintenance. Three rates were calibrated to P6, all confidence C
([parameters.md §8](parameters.md#8-saliva-plaque-and-diet)): growth on
saliva at up to 0.06 per hour, maintenance at 0.0015 per hour, and
respiration of sugar at 0.0035 per hour. Against them, the wear of 0.5 µm per
hour was chosen so that the day returns the plaque within 5% of where it
started, 62.8 µm against 60.

Other checks:

- **Each brushing** takes 42% of the plaque. The record five minutes after
  07:45 finds it 41.9% thinner, having grown back a little since.
- **The plaque grows back** between the brushings, from 39.5 to 84.3 µm, on
  its meals and on saliva.
- **Replay** is the same on 1 and on 4 threads, for both scenes. The day's
  first hour, with breakfast and a brushing, is a slow test.
- **The example** `python examples/plaque_day.py` runs both scenes and checks
  P6, the brushings, the regrowth and P1 itself, in about three minutes.

## Running the suite

`python -m pytest` runs the pull-request suite, which leaves out tests marked
`slow` (long parameter sweeps). `python -m pytest -m slow` runs those, as the
nightly workflow does, and `python -m pytest -m ""` runs everything. As cases
land, `python -m pytest -m "numerical or invariance or regression"` runs the
fast checks and `python -m pytest -m benchmark` runs the biological
benchmarks. A `marse validate` command will run the same suite from an
installed package.
