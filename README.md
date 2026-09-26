# MARSE: Microbial Adaptability Resource Simulation Engine

[![CI](https://github.com/Y3K5/MARSE/actions/workflows/ci.yml/badge.svg)](https://github.com/Y3K5/MARSE/actions/workflows/ci.yml)
[![Privacy](https://github.com/Y3K5/MARSE/actions/workflows/privacy.yml/badge.svg)](https://github.com/Y3K5/MARSE/actions/workflows/privacy.yml)

MARSE is open-source research software for **reproducible, spatial simulation of
microbial populations and biofilms**. It models species, resources and
environmental conditions on a spatial domain, represents adaptation as explicit
state transitions, applies controlled perturbations, and records the full
provenance of every run so that any result can be replayed.

> **Status: Phase 3, microbial dynamics.** MARSE can define, run, save and
> reproduce a well-mixed batch simulation; solve the steady reaction-diffusion
> profile through a biofilm; and couple the two, so growth follows the local
> concentration rather than a well-mixed average. Letting that growth reshape
> the biofilm over time is next; see the [roadmap](docs/roadmap.md).

Run a simulation and then reproduce it exactly from its manifest:

```bash
# a well-mixed batch culture, evolved over time
marse run examples/experiments/two_species_batch.json -o runs/batch
marse replay runs/batch/manifest.json

# a biofilm, solved to the steady profile its oxygen gradient supports
marse run examples/experiments/biofilm_oxygen_profile.json -o runs/biofilm
marse replay runs/biofilm/manifest.json

# a two-dimensional multispecies ecosystem, with a browser viewer
marse ecosystem examples/experiments/two_species_ecosystem.json -o runs/ecosystem
marse replay runs/ecosystem/manifest.json
```

Adding a `biofilm` block to an experiment changes what it *is*: a batch run
evolves a well-mixed culture over time, while a biofilm run holds the biomass
fixed and solves the depth profile, which has no time axis. All three kinds
produce a manifest that replays them. Replay proves a result can be
reproduced, not that it is right: the ecosystem engine has known defects that
are being fixed before anything is built on it
([validation.md](docs/validation.md#the-two-dimensional-ecosystem-engine-is-not-yet-verified)).

The fix starts with the configuration. In the new format every process must
conserve carbon, nitrogen and electrons exactly, and MARSE derives what a
yield leaves open, such as the oxygen used and the carbon dioxide released. A
process that would make matter from nothing is refused as the file loads. A
network with rates runs in a closed, well-mixed box, and every run proves its
own balance:

```bash
marse check examples/networks/glucose_cross_feeding.json
marse run examples/networks/glucose_cross_feeding.json -o runs/network
marse replay runs/network/manifest.json
```

[docs/networks.md](docs/networks.md) describes the format. Transport, and with
it biofilms, comes next.

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
python examples/biofilm_profile.py     # solved oxygen profile through a biofilm
python examples/stratified_growth.py   # why a thick biofilm grows no faster than a thin one
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
  schemas/         configuration schema v2: formulas and balanced reaction networks
  spatial/         domains, grids, diffusion; solute properties
  microbes/        growth kinetics, cardinal models, niches, genotypes, dose responses
  biofilm/         biomass, matrix, maturation
  ecosystem/       the 2-D multispecies engine, frame store and viewer
  adaptation/      state transitions and decision policies
  interventions/   perturbations
  analysis/        calibration, uncertainty and sensitivity, ensembles
  evidence/        culture conditions and measurements, with their sources
  experimental/    outside the v1.0 claims (host/: host-pressure primitives)
  validation/      analytical references, benchmark and regression cases
examples/          runnable reference calculations
  experiments/     experiment configurations for `marse run` and `marse ecosystem`
  networks/        reaction networks for `marse check` and `marse run` (configuration schema v2)
tests/             test suite
docs/              theory, parameters, specification, architecture, validation, roadmap
tools/             repository tooling (privacy and repository guards)
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
- [docs/networks.md](docs/networks.md): reaction networks, the configuration
  format in which every process must conserve carbon, nitrogen and electrons.
- [docs/validation.md](docs/validation.md): benchmark definitions.
- [docs/modeling-landscape.md](docs/modeling-landscape.md): how the field
  models microbial growth, how natural conditions differ from laboratory ones,
  and what both imply for MARSE's design.
- [docs/roadmap.md](docs/roadmap.md): development phases.

## Contributing, privacy and security

- [CONTRIBUTING.md](CONTRIBUTING.md): workflow, tests and scientific standards.
- [GOVERNANCE.md](GOVERNANCE.md): who may merge and release, and the settings that enforce it.
- [docs/repository-setup.md](docs/repository-setup.md): applying those settings, step by step.
- [PRIVACY.md](PRIVACY.md): what never enters this repository, and how that is enforced.
- [SECURITY.md](SECURITY.md): report vulnerabilities or leaked data privately.

## Citation and license

Citation metadata is in [CITATION.cff](CITATION.cff) (GitHub's "Cite this
repository" button reads it). MARSE is licensed under the
[Apache License 2.0](LICENSE).
