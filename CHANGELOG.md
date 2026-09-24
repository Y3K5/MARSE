# Changelog

All notable changes to MARSE are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and versions follow
[Semantic Versioning](https://semver.org/). Any change that alters simulation
results for the same manifest is always called out.

## [Unreleased]

### Added

- **Growth coupled to the gradient (Phase 3, first increment).**
  `biofilm/biomass.py` joins the spatial solute field to the growth kinetics:
  uptake capacity follows the local biomass, the solute field is solved against
  it, and each depth then grows at the rate its own concentration supports.
  Reports per-depth growth rates, the active-zone thickness and the
  biomass-weighted mean growth rate.
- The diffusion solver accepts a per-node uptake capacity, not only a uniform
  one, which is what makes that coupling possible.
- `oxygen_diffusivity_um2_per_h`, because combining a per-second diffusivity
  with per-hour growth rates understates penetration sixtyfold while still
  producing a plausible-looking number. A test pins the size of that error.
- `examples/stratified_growth.py`: identical cells at uniform density, and the
  gradient alone decides which of them grow.

- **Spatial transport (Phase 2).** `spatial/domain.py` provides a uniform
  one-dimensional grid through the depth of a slab, with an imposed surface
  concentration and an impermeable base. `spatial/diffusion.py` solves the
  steady reaction-diffusion equation with Monod uptake by Newton iteration on
  a tridiagonal system, and reports the concentration profile, surface flux,
  total uptake, penetration depth and anoxic fraction.
- Validation case V3 now passes for the one-dimensional steady state, checked
  against three analytical limits that bracket real Monod uptake — no uptake,
  first order and zero order — plus a flux balance and a convergence study
  confirming the discretisation is second order.
- `first_order_profile` and `first_order_surface_flux` analytical references,
  written so that a thick slab cannot overflow `cosh`.
- `examples/biofilm_profile.py`: solves the oxygen profile through a biofilm
  and compares the result with the closed-form penetration depth.

- `.github/rulesets/`: branch and tag protection committed as importable JSON,
  so a change to who may write to the default branch is reviewable like any
  other change. The branch ruleset requires a pull request with all seven CI
  and privacy checks passing and blocks force pushes and deletion; the tag
  ruleset stops a `v*` release tag being moved or deleted.

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
