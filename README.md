# MARSE: Microbial Adaptability Resource Engine

[![CI](https://github.com/Y3K5/MARSE/actions/workflows/ci.yml/badge.svg)](https://github.com/Y3K5/MARSE/actions/workflows/ci.yml)
[![Privacy](https://github.com/Y3K5/MARSE/actions/workflows/privacy.yml/badge.svg)](https://github.com/Y3K5/MARSE/actions/workflows/privacy.yml)

MARSE is open-source research software for **reproducible, spatial simulation of
microbial populations and biofilms**. It models species, resources and
environmental conditions on a spatial domain, represents adaptation as explicit
state transitions, applies controlled perturbations, and records the full
provenance of every run so that any result can be replayed.

> **Status: Phase 0, specification.** This repository holds the design, the
> package layout and the development infrastructure. There is no simulation
> engine yet; see the [roadmap](docs/roadmap.md).

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
  core/            run loop, state, configuration, provenance
  schemas/         validated canonical objects
  spatial/         domains, grids, diffusion and transport
  microbes/        growth, resource use, adhesion, interactions
  biofilm/         biomass, matrix, maturation
  adaptation/      state transitions and decision policies
  interventions/   perturbations
  validation/      benchmark and regression cases
tests/             test suite
docs/              specification, architecture, validation plan, roadmap
tools/             repository tooling (privacy guard)
```

Planned: `examples/` for reproducible reference simulations and `paper/` for the
software paper.

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

## Contributing, privacy and security

- [CONTRIBUTING.md](CONTRIBUTING.md): workflow, tests and scientific standards.
- [PRIVACY.md](PRIVACY.md): what never enters this repository, and how that is enforced.
- [SECURITY.md](SECURITY.md): report vulnerabilities or leaked data privately.

## Citation and license

Citation metadata is in [CITATION.cff](CITATION.cff) (GitHub's "Cite this
repository" button reads it). MARSE is licensed under the
[Apache License 2.0](LICENSE).
