# Changelog

All notable changes to MARSE are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and versions follow
[Semantic Versioning](https://semver.org/). Any change that alters simulation
results for the same manifest is always called out.

## [Unreleased]

### Added

- **Repository hardening.** `tools/repo_guard.py` runs on every commit and in
  CI, and fails if a GitHub Action is not pinned to a commit SHA, a workflow
  uses `pull_request_target`, a workflow takes write permissions that were not
  declared, untrusted text is interpolated into a shell command, a publishing
  step appears without an approval gate, a required file is missing, the
  documented package layout is broken, or unexpected files accumulate at the
  repository root. In CI it runs directly rather than through pre-commit, so
  disabling it requires editing a workflow.
- `GOVERNANCE.md`: who may merge and release, how a release is made, and the
  GitHub settings that enforce access control.
- `.github/CODEOWNERS`: the maintainer's review is requested on every path,
  with the workflow, tooling and policy paths named separately.

- **Simulation kernel (Phase 1).** A well-mixed batch culture of one or more
  organisms competing for a single limiting substrate can now be defined, run,
  saved and reproduced:
  - `marse run <experiment.json>` writes a trajectory and a run manifest;
    `marse replay <manifest.json>` rebuilds the configuration, verifies its
    checksum, re-runs and reports whether the results match.
  - Validated configuration with units in every field name; unknown fields and
    out-of-range values are refused with the offending field named.
  - Deterministic RK4 clock, checkpoints, and a substrate balance checked every
    step. A diverging run or a timestep too large for the configured rates is
    reported rather than silently clamped.
  - Name-keyed random streams, so adding a provider cannot perturb the stream
    any other provider sees.
  - Run manifests recording the configuration, its SHA-256 checksum, the seed,
    versioned model names and software versions — and no username, hostname or
    absolute path.

- `docs/modeling-landscape.md`: survey of established microbial and biofilm
  simulators, how growth in natural environments differs from laboratory
  growth, and the consequences for MARSE's design.
- A growth threshold: `net_growth_rate` and `minimum_substrate_concentration`
  implement `S_min`, below which maintenance consumes all uptake and there is
  no net growth. Bare Monod kinetics predicts growth at any positive substrate
  concentration, which is wrong at environmental concentrations.
- Growth kinetics: Monod, Haldane–Andrews substrate inhibition, Pirt uptake
  with maintenance, Luedeking–Piret product formation, and the Baranyi–Roberts
  growth curve with lag, evaluated in a numerically stable form.
- Secondary models: cardinal temperature (CTMI) and pH (CPM) models and the
  Ratkowsky square-root model, following the gamma concept.
- Oxygen solubility (Benson–Krause) and diffusivity (Han–Bartels) correlations,
  which refuse to extrapolate outside their stated validity ranges.
- Analytical reference solutions for validation cases V1–V5, each checked in
  the test suite against an independent numerical solution.
- `docs/theory.md` (every equation, its assumptions, numerical methods and
  limitations) and `docs/parameters.md` (values with units, sources and
  explicit confidence levels).
- Runnable examples: `batch_growth.py` and `oxygen_penetration.py`.
- Package layout following the planned architecture, and a `marse --version`
  command.
- Scientific specification, architecture, validation plan and roadmap in `docs/`.
- Privacy protection: local git hooks (gitleaks, plus a privacy guard for files,
  commit messages, git identity and outgoing commits) and a CI workflow that
  checks every file and every commit and scans the full history for secrets.
- CI running the pre-commit hooks and the tests on Linux (Python 3.12 to 3.14)
  and Windows.
- Apache-2.0 license and citation metadata.
