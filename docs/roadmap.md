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
   neutralization, and deterministic immune-cell action primitives are
   implemented; receptor kinetics, host tissue, and clinical models remain
   out of scope.)*

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
