# MARSE: Microbial Adaptability Resource Engine

[![CI](https://github.com/Y3K5/MARSE/actions/workflows/ci.yml/badge.svg)](https://github.com/Y3K5/MARSE/actions/workflows/ci.yml)
[![Privacy](https://github.com/Y3K5/MARSE/actions/workflows/privacy.yml/badge.svg)](https://github.com/Y3K5/MARSE/actions/workflows/privacy.yml)

MARSE is open-source research software for **reproducible, spatial simulation of
microbial populations and biofilms**. It models species, resources and
environmental conditions on a spatial domain, represents adaptation as explicit
state transitions, applies controlled perturbations, and records the full
provenance of every run so that any result can be replayed.

> **Status: Phase 1, simulation kernel.** MARSE can define, run, save and
> reproduce a well-mixed batch simulation. Space, diffusion and biofilm
> structure are Phases 2–4; see the [roadmap](docs/roadmap.md).

Run a two-species batch culture and then reproduce it exactly from its manifest:

```bash
marse run examples/experiments/two_species_batch.json -o runs/demo
marse replay runs/demo/manifest.json
```

The manifest records the configuration, its SHA-256 checksum, the random seed,
the versioned models used and the software versions — and deliberately records
no username, hostname or absolute path, so it is safe to attach to a paper or
an issue.

The mathematics underneath is implemented and tested: growth kinetics, cardinal
temperature and pH models, oxygen solubility and diffusivity, and analytical
reference solutions that the validation suite checks against independent
numerical integration. See [docs/theory.md](docs/theory.md) for the equations
and [docs/parameters.md](docs/parameters.md) for sourced parameter values.

```bash
python examples/batch_growth.py        # Monod growth vs its closed-form solution
python examples/oxygen_penetration.py  # why active biofilms are stratified
```

## What MARSE v1.0 will do

- Define microbial species and strains with parameterized growth and resource use.
- Simulate spatial domains with temperature, pH, oxygen and resource fields.
- Model resource diffusion independently of microbial state.
- Represent biomass, adhesion, local competition and cooperation, and simple
  extracellular-matrix effects.
- Run multispecies simulations through explicit interaction rules.
- Model adaptation as configurable state transitions (planktonic, attached,
  stressed, slow-growing, tolerant).
- Apply controlled perturbations and compare them with baseline runs.
- Save every run with its parameters, random seed, versions and assumptions, so
  it can be replayed exactly.
- Ship a validation suite of benchmark and regression cases.

MARSE is a simulation framework. It does not model the human immune system,
predict clinical or vaccine outcomes, or produce experimental evidence; the
[specification](docs/specification.md) sets out these boundaries.

## Design principles

- **Small deterministic core, replaceable scientific providers.** The core owns
  state, time, seeds and provenance; providers own the biology.
- **Explicit scales.** Modules exchange data only through declared state
  variables, and every output records its scale, resolution, rules,
  assumptions, parameter sources and uncertainty.
- **Provenance first.** Every run can be replayed from its manifest.
- **Validation before new biology.** Benchmark cases are defined before the
  models they test.
- **Private by design.** No telemetry, no network access, and run manifests
  never record usernames, hostnames or absolute paths.

## Repository layout

```text
src/marse/         the Python package
  core/            run loop, state, configuration, seeds, provenance
  schemas/         validated canonical objects
  spatial/         domains, grids, diffusion; solute properties
  microbes/        growth kinetics, cardinal models, interactions
  biofilm/         biomass, matrix, maturation
  adaptation/      state transitions and decision policies
  interventions/   perturbations
  validation/      analytical references, benchmark and regression cases
examples/          runnable reference calculations
  experiments/     experiment configurations for `marse run`
tests/             test suite
docs/              theory, parameters, specification, architecture, validation, roadmap
tools/             repository tooling (privacy guard)
```

Planned: `paper/` for the software paper.

## Development setup

Requires Python 3.12 or newer.

```bash
git clone https://github.com/Y3K5/MARSE.git
cd MARSE
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -e ".[dev]"
pre-commit install
python -m pytest
marse --version
```

Before your first commit, complete the one-time setup in
[PRIVACY.md](PRIVACY.md#one-time-setup).

## Documentation

- [docs/theory.md](docs/theory.md): every equation, its assumptions, and the
  numerical methods used.
- [docs/parameters.md](docs/parameters.md): parameter values with units,
  sources and explicit confidence levels.
- [docs/specification.md](docs/specification.md): scientific scope of v1.0,
  and what it does not claim.
- [docs/architecture.md](docs/architecture.md): core abstractions and
  extension contracts.
- [docs/validation.md](docs/validation.md): benchmark definitions.
- [docs/modeling-landscape.md](docs/modeling-landscape.md): how the field
  models microbial growth, how natural conditions differ from laboratory ones,
  and what both imply for MARSE's design.
- [docs/roadmap.md](docs/roadmap.md): development phases.

## Contributing, privacy and security

- [CONTRIBUTING.md](CONTRIBUTING.md): workflow, tests and scientific standards.
- [PRIVACY.md](PRIVACY.md): what never enters this repository, and how that is enforced.
- [SECURITY.md](SECURITY.md): report vulnerabilities or leaked data privately.

## Citation and license

Citation metadata is in [CITATION.cff](CITATION.cff) (GitHub's "Cite this
repository" button reads it). MARSE is licensed under the
[Apache License 2.0](LICENSE).
