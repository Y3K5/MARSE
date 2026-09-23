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
| V1 | Single-species unrestricted growth | The chosen growth law and carrying-capacity behaviour | 3 | Planned |
| V2 | Resource-limited growth | Expected saturation or starvation behaviour | 3 | Planned |
| V3 | Diffusion only | The numerical diffusion method, separately from microbial rules | 2 | Planned |
| V4 | Attachment and biofilm initiation | Transition from planktonic or seeded biomass to attached growth | 4 | Planned |
| V5 | Two-species competition | Expected dominance and coexistence regimes | 3 | Planned |
| V6 | Cross-feeding | Explicit beneficial exchange shifts the equilibrium as designed | 3 | Planned |
| V7 | Environmental perturbation | Recovery or adaptation after a resource, pH or oxygen shift | 4 | Planned |
| V8 | Reproducibility | Re-running a saved manifest reproduces outputs within stated tolerances | 5 | Planned |

## Running the suite

As cases land, `python -m pytest -m "numerical or invariance or regression"`
runs the fast checks and `python -m pytest -m benchmark` runs the biological
benchmarks. A `marse validate` command will run the same suite from an
installed package.
