# Roadmap

MARSE is developed in the open, in phases. Each phase has an exit criterion,
and the next phase starts once it is met. Progress is tracked in issues and
milestones.

**Current phase: 3, microbial dynamics.** A well-mixed batch simulation can be
defined, run, saved and reproduced exactly from its manifest (Phase 1); the
one-dimensional reaction-diffusion solver is validated against three analytical
limits (Phase 2, validation case V3); and growth is now coupled to the solute
gradient, so each depth grows at the rate its own local concentration supports.
What remains in Phase 3 is letting that growth change the biomass distribution
over time, which is where the biomass-spreading choice of
[modeling-landscape.md](modeling-landscape.md#2-the-biomass-spreading-decision)
has to be made.

Alongside that core, a two-dimensional multispecies engine (`marse.ecosystem`)
now runs, with a browser viewer and a periodontal example study. It is not yet
verified, and a review found defects in its mass balance, transport and
reproducibility ([validation.md](validation.md#the-two-dimensional-ecosystem-engine-is-not-yet-verified)).
Fixing them comes before anything new; see
[the order of work](#order-of-work-correctness-first).

| Phase | Deliverable | Exit criterion |
|---|---|---|
| 0. Specification | Schemas, model assumptions, benchmark plan | Core objects and scientific scope settled enough to implement |
| 1. Simulation kernel | Time and state engine, configuration, checkpoints | Deterministic replay and unit tests pass |
| 2. Space and resources | Domain, neighbourhoods, diffusion and transport | Analytical and numerical benchmarks pass |
| 3. Microbial dynamics | Growth, resource use, adhesion, interactions | Single- and two-species cases are reproducible |
| 4. Biofilm and adaptation | Matrix and maturation; explicit state transitions | Perturbation benchmarks and the sensitivity suite pass |
| 5. Reproducibility | Manifests, provenance, exports, replay | An independent example can be reproduced from the repository |
| 6. Reference analyses | Three to five documented research-style cases | Results and limitations are documented |
| 7. Decision-policy adapters (optional) | Bounded `DecisionPolicy` plugins | No core dependency; the deterministic baseline is kept |
| 8. v1.0 release | Documentation, tests, license, stable API | Feature complete |
| 9. Software paper | Short paper and archived release | Public development history; ready for review |

## Build order

1. Repository with license, contribution guide, issue templates and a
   scientific specification. *(Done.)*
2. Validated configuration models for `Organism`, `Environment`,
   `ResourceField`, `PopulationState`, `Interaction`, `Perturbation` and
   `SimulationRun`. *(Done for the non-spatial subset: organism, environment,
   substrate and run. `ResourceField`, `Interaction` and `Perturbation` arrive
   with Phases 2–4.)*
3. Deterministic simulation clock, state checkpoints, random-seed handling and
   the run manifest. *(Done.)*
4. A simple spatial domain and diffusion solver, validated on its own.
   *(Done: `spatial/domain.py` and `spatial/diffusion.py`, validated against
   no-uptake, first-order and zero-order analytical limits, a flux balance and
   a second-order convergence study.)*
5. One microbial growth model and one resource-consumption model.
6. Local interaction and adhesion rules sufficient for a first biofilm
   demonstration.
7. Three to five benchmark cases, before any more biology.
8. Visualization and export, once the state and validation layers are stable.
9. The optional `DecisionPolicy` interface, with `RuleBasedPolicy` as the
   reference.
10. Experimental external policy adapters, compared against the deterministic
    policies on bounded decisions.

## Order of work: correctness first

The ecosystem engine grew faster than its verification. Until the defects in
[validation.md](validation.md#the-two-dimensional-ecosystem-engine-is-not-yet-verified)
are fixed, **no new mechanism is added**: work is limited to fixes,
verification and documentation. After that, a mechanism lands only with its
equations in `theory.md`, units in its configuration field names, a
conservation or invariance test, a verification or validation case that runs
the simulator, graded parameters in `parameters.md`, and a `CHANGELOG` entry.

| Stage | Work | Done when |
|---|---|---|
| 0. Stabilise ✓ | Honest docs and tests; each known defect written as an expected-failure test; long sweeps moved to a nightly run | Pull-request CI is fast and every defect has a test |
| 1. Provenance and infrastructure | Ecosystem runs go through the core manifest and replay *(done)*; frames recorded on demand and streamed to disk instead of every frame held in memory *(done)*; the package layout (`ecosystem/`, `evidence/`, `analysis/`, `experimental/host/`). The configuration is redesigned in Stage 2, not here. | `marse replay` reproduces an ecosystem run exactly *(met)*; memory no longer grows with run length *(met)* |
| 2. The material core | Conservation by construction, as in the IWA biofilm models. Components have a composition (C, N and electrons). Processes are rows of a stoichiometric (Gujer) matrix, and a continuity check refuses at load time any process that creates matter. One rate per process drives every component it touches. Reactions use a conservative, positive (Patankar-type) integrator; transport is conservative implicit finite-volume diffusion at physical diffusivities. Each species has an oxygen role, variants are heritable lineages, and space is shared. A ledger checks the balance every step. Configuration schema v2 replaces v1, which is removed, and the examples and periodontal study are rebuilt on it. Built in increments: (2a) the v2 reaction-network language, whose processes are checked for continuity when they load *(done: [networks.md](networks.md))*; (2b) a well-mixed engine with the per-step ledger *(done)*; (2c) transport; (2d) oxygen roles and lineages; (2e) the examples and study rebuilt, v1 removed. | Every known-defect test passes; a closed box conserves C, N and electrons to 1e-12 *(met in 2b)* |
| 3. Compartments and elements | Dead biomass, lysis and hydrolysis, extracellular matrix and inert pools; porosity and effective diffusivity from the matrix; detachment recorded as an export; named decay products; sulphur and phosphorus | Open-system benchmarks: the import at steady state matches the analytical surface flux |
| 4. Verification | 2-D transport against `point_source_diffusion_2d`; the 2-D solver reducing to V3; manufactured solutions; permutation invariance; the R\* rule through the simulator | Convergence at design order; cross-platform replay in CI |
| 5. External validation | IWA BM1 (V9) and BM3 (V10); oxygen microprofiles | Published reference solutions reproduced |
| 6. Flagship study | The periodontal biofilm rebuilt with graded parameters, named metabolite exchanges and a published in-vitro calibration target | A falsifiable, evidence-linked prediction |
| 7. Scale | Year-long runs and large ensembles, streamed and resumable | Annual runs on a laptop |
| 8. Publication | JOSS paper on the verified software; a validation paper on Stages 5–6 | Submitted |

The long-term plan below is sequenced after these stages, not beside them.

## Long-term scientific platform plan

The post-Phase-3 plan is to make every result evidence-linked and falsifiable:

1. **Evidence-to-experiment compiler** — compile only context-compatible
   culture records into validated configurations; preserve source IDs,
   assumptions, unresolved values, and dataset checksums. *(Broth slice
   implemented; spatial formats remain explicit.)*
2. **Provider-based ecosystem kernel** — separate transport, reactions,
   capability responses, biomass, events, mechanics, and diagnostics behind
   versioned contracts with fixed update order.
3. **Calibration and inverse modelling** — fit growth, CFU, colony-area,
   inhibition, nutrient, and oxygen observations with replicate-aware
   uncertainty and identifiability checks. *(Growth-curve screening fits
   implemented; confidence intervals, censoring, and full inverse modelling
   remain.)*
4. **Uncertainty and sensitivity** — propagate measurement, parameter,
   numerical, and model-form uncertainty through reproducible ensembles.
   *(Bounded sampling, summary quantiles, and rank-based screening implemented;
   posterior inference and model-form comparison remain.)*
5. **Benchmark organisms and experiments** — maintain a small source-linked
   validation panel before expanding mechanism count. *(Analytical benchmark
   registry and comparison contracts implemented; curated organism cases
   remain.)*
6. **Adaptive spatial behaviour** — add chemotaxis, quorum signals, adhesion,
   detachment, morphology, and phenotype switching after the provider contracts
   are stable. *(Optional conservative chemotaxis, explicit diffusing quorum
   signals, quorum-triggered phenotype switching, and opt-in surface
   adhesion/detachment are implemented; morphology remains.)*
7. **Genotype-to-capability maps** — connect genotype to regulated protein
   abundance and traits only where data constrain the mapping. *(Explicit
   allele-to-capability parameter modifiers are implemented; sequence
   interpretation and regulated abundance remain future work.)*
8. **Scientific observability** — expose conservation residuals, validity
   limits, evidence compatibility, convergence, and provenance for every run.
9. **Immune and molecular extensions** — add explicit effector pressure and
   molecular neutralization primitives. *(Field-level deterministic pressure
   neutralization, deterministic immune-cell action primitives, and optional
   ecosystem pressure integration are implemented; receptor kinetics, host
   tissue, and clinical models remain out of scope.)*
10. **Action and resource economy** — make typed actions consume explicit
    resources and require declared capabilities. *(Reusable budget and
    prerequisite accounting and immune-agent budget enforcement are
    implemented; spatial resource coupling remains.)*

## First publishable milestone

A researcher can install MARSE, run a documented single- and multispecies
biofilm simulation, perturb the environment, reproduce the run from its
manifest, inspect parameter provenance, and run the validation suite.

## Publication

The software paper will be short and focused on the software: summary,
statement of need, design, validation, research use and limitations, and
availability. Detailed theory stays in the documentation. The target venue is
the Journal of Open Source Software, after v1.0 and at least six months of
public development. Journal requirements will be rechecked before submission.
