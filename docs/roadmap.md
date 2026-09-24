# Roadmap

MARSE is developed in the open, in phases. Each phase has an exit criterion,
and the next phase starts once it is met. Progress is tracked in issues and
milestones.

**Current phase: 1, simulation kernel.** A well-mixed batch simulation can be
defined, run, saved and reproduced exactly from its manifest. Space and
diffusion are next.

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
5. One microbial growth model and one resource-consumption model.
6. Local interaction and adhesion rules sufficient for a first biofilm
   demonstration.
7. Three to five benchmark cases, before any more biology.
8. Visualization and export, once the state and validation layers are stable.
9. The optional `DecisionPolicy` interface, with `RuleBasedPolicy` as the
   reference.
10. Experimental external policy adapters, compared against the deterministic
    policies on bounded decisions.

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
