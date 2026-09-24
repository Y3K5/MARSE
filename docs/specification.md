# Scientific specification

*Phase 0 draft. This document fixes the scientific scope of MARSE v1.0 and
changes only through reviewed pull requests.*

## Purpose

MARSE is research software for reproducible, spatial simulation of microbial
populations and biofilms under defined environmental conditions. It is built
for questions such as:

- How does resource limitation shape biofilm growth and structure?
- Under which interaction parameters do two species coexist, and when does one
  dominate?
- How does a community respond to, and recover from, a nutrient, pH or oxygen
  shift?
- Which parameters drive a result, and how robust is a conclusion to
  parameter uncertainty?

Every answer MARSE gives is conditional on the simulated conditions and the
stated model assumptions.

## Scope of v1.0

MARSE v1.0 consists of a simulation kernel; species, environment and resource
models; spatial biofilm dynamics; configurable interactions; perturbations;
reproducibility; and a validation suite.

At v1.0 a researcher can install MARSE, run a documented single- and
multispecies biofilm simulation, perturb the environment, reproduce the run
from its manifest, inspect where every parameter came from, and run the
validation suite.

### Core capabilities

1. Species and strains with parameterized growth and resource use.
2. Spatial domains with boundary conditions; temperature, pH, oxygen and
   resource fields; a defined time resolution.
3. Resource diffusion or transport, updated independently of microbial state.
4. Biomass, adhesion, local competition and cooperation, and simple
   extracellular-matrix effects.
5. Multispecies simulations through explicit interaction rules, never hidden
   heuristics.
6. Adaptation as configurable state transitions (for example planktonic,
   attached, stressed, slow-growing or tolerant) rather than genomic evolution.
7. Controlled perturbations, compared against baseline trajectories.
8. A complete record of every run: parameters, random seed, software and model
   versions, assumptions and output summaries.
9. Benchmark and regression cases that demonstrate numerical stability and
   selected biological fidelity.

### Deferred beyond v1.0

- Detailed simulation of host immune cells.
- Atomistic molecular dynamics as an integrated requirement.
- Vaccine efficacy prediction, and de novo antigen or molecular design.
- Clinical decision support.
- Learned adaptive policies that change biological rules without an audit trail.
- Claims of organism-wide or whole-host "digital twin" fidelity.

The architecture leaves interfaces for richer metabolism, host interaction,
molecular providers and multiscale coupling
([architecture](architecture.md#extension-points-after-v10)), but v1.0 does
not implement them.

## Claims

| MARSE v1.0 may claim | MARSE v1.0 must not claim |
|---|---|
| A reproducible framework for spatial microbial and biofilm simulation | A complete simulated human immune system |
| Configurable resource, environmental and interaction models | Clinically predictive vaccine efficacy |
| Multispecies and perturbation experiments under explicit assumptions | Whole-body physiological realism |
| Agreement with selected published or experimental behaviours, as benchmarked | That molecular dynamics alone explains colony-scale behaviour |
| An extensible architecture for future multiscale coupling | That model- or AI-generated output is experimental evidence |

## Scales

| Scale | Method in v1.0 | Later extension |
|---|---|---|
| Molecular | Not included | Docking and molecular dynamics provider interface |
| Cell or microbial unit | Agent or local biomass state | Richer phenotype and metabolic models |
| Biofilm | Spatial occupancy, matrix and diffusion | Mechanical growth and morphology models |
| Community | Explicit multispecies interactions | Genome-scale metabolism and ecology |
| Host | Deferred | Simplified host-cell and immune-state model |
| Intervention | Abstract perturbation pressure | Molecularly grounded interventions |

## Canonical objects

| Object | Responsibility |
|---|---|
| `Organism` | Species or strain identity, with its growth, resource and phenotype parameters |
| `Environment` | Physical and chemical context, and boundary conditions |
| `ResourceField` | Spatially varying nutrients, oxygen, metabolites or other modelled quantities |
| `PopulationState` | Biomass or cell count, position, phenotype and local state |
| `Interaction` | An explicit rule or function connecting organisms or state variables |
| `Perturbation` | A controlled change applied at a specified time or condition |
| `SimulationRun` | Configuration, seed, versions, checkpoints, outputs and provenance |
| `ValidationCase` | A known behaviour, target metric, tolerance and source |

The equations that implement these objects, with their assumptions and
limitations, are in [theory.md](theory.md); their numerical values and sources
are in [parameters.md](parameters.md).

## Modelling rules

1. **Explicit scales.** Modules connect only through declared state variables
   and interfaces. Mismatched scales are never hidden behind one composite
   "realism" score.
2. **Self-describing outputs.** Every output records its biological scale,
   spatial and temporal resolution, equations or rule set, assumptions,
   parameter sources and uncertainty.
3. **Sourced parameters.** Every parameter carries units, a source (publication
   or dataset) and a confidence note. Physical constants and biological rates
   come from data and model definitions, never from a learned policy.
4. **Controlled randomness.** All stochasticity comes from seeded random
   streams owned by the core. The same manifest reproduces the same run, or,
   for stochastic ensembles, the same results within stated tolerances.
5. **Simulation is not evidence.** Conclusions are stated in terms of the
   simulated conditions and assumptions, not as clinical or experimental fact.

## Open questions for Phase 0

To be settled in issues before Phase 1 starts:

- Spatial representation for v1.0: a lattice of biomass densities, individual
  agents, or a hybrid? **Partly answered** by
  [modeling-landscape.md §2](modeling-landscape.md#2-the-biomass-spreading-decision):
  the choice of biomass-spreading mechanism measurably changes conclusions
  about competition and cooperation, so it should be a declared, swappable
  provider rather than a single fixed choice.
- Dimensionality: two dimensions first, three later?
- Units: one internal system (for example µm, minutes and mmol) with
  conversion at the configuration boundary?
- Diffusion solver: explicit finite differences with a stability check, or an
  implicit scheme from the start?
- Configuration format, and versioning of the experiment and manifest schemas.
- Output formats (for example Parquet trajectories with JSON or YAML
  summaries) and the dependencies they bring.
