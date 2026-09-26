# Changelog

All notable changes to MARSE are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and versions follow
[Semantic Versioning](https://semver.org/). Any change that alters simulation
results for the same manifest is always called out.

## [Unreleased]

### Changed

- **The package follows its documented layout.** Modules that had been added
  at the package root now live in subpackages:
  - `niche`, `genotype` and `additives` are in `marse.microbes`;
  - `science` is in `marse.evidence`;
  - `calibration`, `uncertainty` and `ensemble` are in `marse.analysis`;
  - `immune` and `actions` are in `marse.experimental.host`, which makes the
    specification's "experimental, outside the v1.0 claims" visible in every
    import.

  `tools/repo_guard.py` now requires the new subpackages, and the module maps
  in the README and `docs/architecture.md` match the tree.

- **Frames are recorded on demand and streamed to disk.** `marse ecosystem`
  writes its frames to a `frames/` store: one NumPy `.npy` file per field,
  in single precision, plus an `index.json`. This replaces `frames.json`, and
  memory stays flat however long the run. It stores about 200 evenly spaced
  frames by default; `--frame-interval-h` sets the spacing. On the two-species
  example the outputs fall from 314 MB to 43 MB, the run from 25 s to 3 s, and
  peak memory from 952 MB to 235 MB. **Breaking:** `marse ecosystem` and
  `write_viewer` no longer write `frames.json`. `EcosystemResult.write_frames`
  remains for library use until configuration schema v2 replaces this engine.
- `run()` takes `frame_every` (record every so many steps plus the last;
  `None` for only the first and last) and `sink` (receive each frame as it is
  made). Recording only observes: a test requires a bit-identical final state
  for every setting. The niche maps a frame shows are computed only for frames
  that are recorded. They had been computed every step, which was 37% of a
  periodontal run's time.
- The ensemble and both periodontal sweeps read only the final state, so they
  now record no intermediate frames.
- The viewer is self-contained. It embeds up to 100 evenly spaced frames,
  always the first and last, with fields rounded to four significant digits
  and totals from the unrounded values.

- Manifest format version 2 adds `kind` (`batch`, `biofilm_profile` or
  `ecosystem`), and replay rebuilds the configuration with the matching parser,
  refusing a manifest whose configuration contradicts its kind. Version 1
  manifests remain readable. Batch and biofilm results, and their `run_id`s,
  are unchanged.
- The roadmap's stages follow revision 2 of the plan: Stage 2 builds the
  conserving material core (stoichiometric processes with a load-time
  continuity check and an every-step ledger) and replaces configuration
  schema v1.

- **The two-dimensional ecosystem engine is documented as unverified.** A
  review that ran it found nine defects. Among them: uptake follows potential
  rather than actual growth, production has no source, the examples give oxygen
  a diffusivity about 4×10⁵ below its physical value, and the periodontal
  anaerobes are modelled as needing oxygen. `docs/validation.md` lists them with
  measurements. Each is a test in `tests/test_ecosystem_known_defects.py` that
  asserts the correct behaviour and is marked as an expected failure, so the fix
  is what flips it. No simulation result changes in this release.
- `docs/roadmap.md` puts correctness first: no new mechanism until the known
  defects are fixed, and after that a mechanism lands only with its equations,
  units, conservation test, verification case, graded parameters and changelog
  entry.
- `docs/specification.md` marks the immune and action primitives as
  experimental and outside the v1.0 claims, and names the periodontal biofilm
  as the flagship application, bounded as a research example.
- `docs/architecture.md` no longer calls the two-dimensional transport operator
  validated. It has not been checked against an analytical solution.
- The benchmark registry's first test scored each case against its own expected
  values, so it could not fail. It is now named for what it does check (the
  comparison machinery). A new test scores the first-order oxygen case against
  the validated one-dimensional solver and confirms second-order convergence.
- Long parameter sweeps are marked `slow`, left out of pull-request CI and run
  nightly (`.github/workflows/nightly.yml`). The pull-request suite drops from
  about nine minutes to about one.
- The project is named "Microbial Adaptability Resource Simulation Engine"
  throughout (citation metadata, package metadata and CLI), matching the README
  and the acronym.

### Deprecated

- The old import paths (`marse.niche`, `marse.genotype`, `marse.additives`,
  `marse.science`, `marse.calibration`, `marse.uncertainty`, `marse.ensemble`,
  `marse.immune`, `marse.actions`) still work for one release. Each emits a
  `DeprecationWarning` and exports the same objects as its new location. A
  test holds every alias to that, and another imports all of MARSE with
  warnings as errors, so nothing inside MARSE uses an alias.

### Added

- **Reaction networks: configuration schema version 2, part one**
  (`marse.schemas`, `marse check`, [docs/networks.md](docs/networks.md)). A
  network lists components, each with a chemical formula, and processes, each
  a row of a stoichiometric matrix. Every process is proven, as it loads, to
  conserve carbon, nitrogen and electrons exactly. The arithmetic is rational,
  on the decimals as written, so there is no tolerance for a leak to hide in.
  A process that would create or destroy matter is refused with an error
  naming the process and the quantity, so known defect 2's "production from
  nothing" cannot be configured in version 2.
  - A growth process states its yield. The coefficients listed in
    `balanced_by` are solved from the balances: typically the electron
    acceptor, carbon dioxide and the nitrogen source, or a fermentation
    product. Water and protons close the oxygen, hydrogen and charge balances
    implicitly.
  - Every key is checked. An unknown key is refused with the likely intended
    one (`yeild_mol_per_mol` → `yield_mol_per_mol`), and a key in the wrong
    unit with the right name (`yield_g_per_g` → `yield_mol_per_mol`).
    Numeric fields carry their unit in their name, and a test enforces it.
  - `marse check NETWORK.json` prints each component's composition and each
    process as a balanced equation. `examples/networks/glucose_cross_feeding.json`
    shows respiration, fermentation and lactate cross-feeding.
  - Verified against textbook degrees of reduction and COD factors, textbook
    reactions, the half-reaction method of Rittmann and McCarty (2001), and
    random networks recounted atom by atom (docs/validation.md, "Stoichiometric
    continuity"). The equations are in docs/theory.md §3.6.

- **Networks run: the version 2 well-mixed engine** (`marse run` on a version
  2 file, `marse.core.well_mixed`). A process with a `rate` runs at
  k × c[proportional_to] × its Monod, inhibition or Haldane factors. One rate
  drives every component the process touches, so uptake is growth divided by
  yield, exactly, and every process acts on the same state at once. These
  are the requirements behind known defects 2 and 7. A second factor for one
  component is refused (defect 4). So is a process consuming what its rate
  does not depend on, unless the file declares the component
  `assumed_in_excess`.
  - Integration is Heun's method with per-process positivity limiting: exactly
    conservative, never negative, with no clipping (the reaction half of
    defect 5). Only the processes consuming a depleted species slow down.
    Substeps adapt to `relative_tolerance` and
    `absolute_tolerance_mol_per_m3`, deterministically.
  - A ledger checks carbon, nitrogen and electrons after every step. A drift
    beyond 1e-9 of the total stops the run with a `ConservationError`, and
    the largest drift is recorded in the manifest (kind `well_mixed`). Runs
    replay bit for bit.
  - Verified against the analytical Monod batch solution (V2) at second order.
    Also: ten thousand steps within 1e-12, a planted leak that is caught,
    3,000 random cyclic networks at steps up to 10⁶ h that stay positive and
    conserve to 5e-16, and order independence
    (docs/validation.md, "The version 2 network engine").
  - The example network now runs: oxygen runs out, and the fermenter turns
    the remaining glucose into lactate.

- **Ecosystem runs are reproducible.** `marse ecosystem` writes a
  `manifest.json` beside its frames and viewer, and `marse replay` reproduces
  the run exactly, as it already did for batch and biofilm runs. The manifest
  records every effective configuration value (defaults included), the seed,
  the provider versions and the engine, plus readable totals and a SHA-256
  digest of the entire final state, so a replay compares every value rather
  than only the totals. The engine is named `ecosystem_v1_unverified` in the
  manifest, because its known defects (docs/validation.md) travel with any
  result it produces. Known defect 9 is fixed, and its test is now an ordinary
  regression test.

- **Two-dimensional multispecies ecosystem engine** (`marse.ecosystem`,
  `marse ecosystem`). It provides nutrient, condition and additive fields with
  explicit transport and fixed-value boundaries; seeded colonies and colony
  spreading; competition coefficients; production into nutrient fields;
  mutation flags; versioned provider contracts; and a self-contained
  interactive HTML viewer with rendering modes and agent-level overlays. It is
  not yet verified: see Changed above.
- Niche capabilities and condition scans (`marse.niche`, `marse niche-scan`),
  coupled to ecosystem growth and shown as viewer layers.
- Chemotaxis up a field gradient; diffusing quorum signals with
  quorum-triggered phenotype switching and hysteresis; surface adhesion and
  detachment.
- Continuous small-molecule additive fields with dose-response effects
  (`marse.additives`).
- Reproducible ecosystem ensembles over parameter grids (`marse.ensemble`,
  `marse ecosystem-batch`).
- An evidence layer for culture conditions and measurements (`marse.science`),
  and a compiler from evidence records to validated configurations.
- Replicate-aware growth-curve calibration (`marse.calibration`), and
  uncertainty propagation with rank-based sensitivity screening
  (`marse.uncertainty`).
- A benchmark registry with comparison metrics
  (`marse.validation.benchmarks`).
- Explicit genotype-to-capability parameter mappings (`marse.genotype`).
- Experimental host-pressure primitives: effector pressure and molecular
  neutralisation, immune-cell action rules and action budgets (`marse.immune`,
  `marse.actions`), with an integrated resistance scenario. These are outside
  the v1.0 claims.
- A three-species periodontal biofilm study with control and environment
  variants, community-control workflows, and a parameter-provenance sidecar
  that marks every value as an uncalibrated placeholder (confidence C).

- **The spatial model is reachable from the kernel.** Until now `marse run`
  could only do well-mixed batch culture: the biofilm code was library-only,
  so MARSE's reproducibility claim did not cover its most interesting output.
  An experiment may now carry a `biofilm` block (thickness, nodes,
  diffusivity), and `marse run` solves the steady depth profile, writes
  `profile.csv` in place of `trajectory.csv`, and records the transport and
  solver versions in the manifest. `marse replay` reproduces it exactly.
- Time settings are refused rather than ignored on a biofilm experiment,
  because a fixed-biomass profile is solved to steady state and has no clock.
- `examples/experiments/biofilm_oxygen_profile.json`.

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
