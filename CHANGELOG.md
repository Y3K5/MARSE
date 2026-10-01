# Changelog

All notable changes to MARSE are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and versions follow
[Semantic Versioning](https://semver.org/). Any change that alters simulation
results for the same manifest is always called out.

## [Unreleased]

### Changed

- **Columns are solved directly.** A one-dimensional run's linear systems are
  now solved exactly, by block-tridiagonal elimination, instead of by
  multigrid. The 1-D example runs three times faster, in the same 348 steps.
  **Results changed**, in their last digits only: the example's totals agree
  with before to about 1e-14. Boxes in two and three dimensions are
  unchanged, bit for bit.
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

- **The roadmap gains the environments programme, and the stages after it.**
  - **Environments (E1–E4):** binding to surfaces (E1, this release); binding
    from surface physics, with natural waters and rock and mineral surfaces
    (E2); soil as a full 3-D pore structure (E3); and geographic scenes (E4).
  - **Saliva, diet and caries** comes next: high- and low-sugar diets
    compared over months to years.
  - **Host biology** (immune, structural and functional cells, and their
    signalling) comes after v1.0, as the specification plans.
- `marse check` says when nothing diffuses, instead of printing an infinite
  explicit step.

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

- **Networks run in space: the version 2 spatial engine** (`marse run` on a
  version 2 file with a `domain`; `marse.core.reactive_transport`).
  - A box of voxels over a surface, in one, two or three dimensions, with the
    bulk liquid held above it.
  - Colonies are placed as hemispheres, either where stated or scattered from
    the run's seed.
  - Dissolved components diffuse at physical diffusivities, the requirement
    behind known defect 1, and every process runs in every voxel. Biomass grows
    in place until Stage 2d lets colonies spread.
  - **Integration** (`marse.core.implicit`) is implicit and coupled, with no
    splitting and no quasi-steady assumption. It uses the two-stage L-stable
    SDIRK method of Alexander (1977): long steps land on the quasi-steady
    profile, and short steps follow transients.
    - The new state is rebuilt from face transfers and process extents, so
      every balance holds to rounding whatever the solver tolerance.
    - A limiter that counts same-step production keeps concentrations
      non-negative without clipping. The update sums what arrives in each
      voxel and what leaves it separately, so a voxel that loses nothing can
      only gain, even in rounding. Traces below the smallest normal number,
      where no relative margin survives rounding, stop giving instead of being
      scaled, so the guarantee holds in floating point.
    - Newton's method is projected, and stays away from the kink of the rate
      laws at zero. Its Jacobian (`rate_jacobian`, analytic) is the slope of the
      rates as evaluated, so a negative value has none.
    - Errors are measured against each component's largest value in the box,
      and a run's first step is estimated from the rates (Hairer, Nørsett and
      Wanner 1993).
  - The ledger books what crosses the top face, computed from the face
    transfers, and checks the box against it after every step. The largest
    residual and the imports are in the manifest (kind `reactive_transport`).
    Runs replay bit for bit.
  - `marse run` writes `totals.csv` (per m² of surface), `frames/`, and `vtk/`
    for ParaView. `marse check` prints the grid, the memory it needs, the
    explicit step it avoids, the multigrid depth, and the diffusion and growth
    time scales.
  - `relative_tolerance` defaults to 1e-4 in space, where the measured error
    is about 5e-6 of each component's peak.
  - Prototyped first in one dimension against criteria set in advance
    (docs/validation.md, "Reactions and transport in space").
  - Measured against its pre-registered performance target, 24 h on a
    64 × 64 × 32 box in under an hour on one core: it took 2 h 17 min
    (docs/validation.md, "Performance"). Most of the excess is the start-up
    transient and the cost per voxel, both taken up in Stage 7.
  - Examples: `examples/networks/surface_biofilm_3d.json`, which runs in about
    five minutes on one core, and a 1-D twin that runs in three seconds.

- **Transport in one, two and three dimensions** (`marse.spatial`), the kernel
  of the spatial engine.
  - `Grid` is a box of cubic voxels over a flat substratum, with height as its
    last axis; one code path serves 1-, 2- and 3-D.
  - `Diffusion` is finite-volume diffusion. Each face flux is computed once and
    applied to both voxels, so matter can only enter or leave through the top
    face, where the bulk liquid is held. The lateral faces are periodic and the
    substratum is impermeable.
  - `ImplicitSystem` is a geometric multigrid solver for the systems an
    implicit step needs, preconditioning GMRES. Each cycle reduces the error
    by a factor of about 0.1 on every grid size tested (8³ to 64×64×32).
  - `write_vti` and `write_pvd` write frames that ParaView opens as a time
    series.
  - The frame store moved to `marse.core.framestore`, so every engine can use
    it; `marse.ecosystem.framestore` still exports it.
  - Verified:
    - a cosine mode is an exact eigenvector of the operator in 1-, 2- and 3-D,
      to rounding;
    - the operator is second order on a 3-D point release, against the new
      reference `point_source_diffusion_3d`;
    - matter is conserved to rounding;
    - a laterally uniform 3-D field diffuses exactly as one column;
    - the multigrid solver agrees with a direct solve.

    See docs/validation.md, "Transport in one, two and three dimensions".

- **Cells bind to surfaces: environments, increment E1**
  ([docs/environments.md](docs/environments.md)). The new modules are
  `marse.spatial.colloids`, `marse.spatial.surface` and
  `marse.microbes.adhesion`. A domain in space can now say what its
  substratum is made of and which cells in the liquid bind to it, so the
  surface is colonized as the run goes instead of being seeded at the start.
  The model is in docs/theory.md §6.4.
  - **Delivery.** Cells diffuse as Brownian spheres (Stokes–Einstein). They
    reach the substratum by the Lévêque flux of a shear flow, the
    Smoluchowski–Levich approximation of flow-chamber studies. A flow
    chamber states its shear. A thin film with a free top, such as saliva on
    a tooth, shears at three times its mean velocity over its thickness.
  - **Binding.** A fraction of the delivered cells binds: the attachment
    efficiency of that species on that material. Bound cells block the area
    around them, and binding stops at the jamming limit of random sequential
    adsorption, a coverage of 0.547 (Feder 1980).
  - **Reversible, then locked.** Bound cells detach, or lock into the
    species' biomass and grow there.
  - **Conservation.** Deposition and detachment cross the substratum face,
    and the ledger books them as imports, beside the top face. Locking
    converts between two components the schema requires to have the same
    formula, as an extra row of the stoichiometric matrix. The exchange has
    an analytic Jacobian, and the positivity limiter scales it like any other
    transfer.
  - **Schema.**
    - Five domain fields go together: `substratum`, `liquid`, `flow`,
      `suspension` and `adhesion`.
    - Patches of material cover the substratum exactly once.
    - Impossible scenes are refused with the reason.
    - A domain without these fields writes exactly what it wrote before, so
      its run id and checksums are unchanged.
  - **Scenes.**
    - `examples/environments/lab/flow_chamber.json`: *S. oralis* on bare and
      saliva-coated glass, for 4 h.
    - `examples/environments/dental/dental_surfaces.json`: enamel, titanium,
      zirconia and acrylic under a salivary film, with *S. oralis* and
      *S. sanguinis*, for 24 h.
    - `examples/surface_adhesion.py`: the four dental materials, hour by
      hour.
  - **Command line.** `marse check` previews delivery and binding on every
    material, and `marse run` writes `surface.csv`.
  - **Verified** against:
    - the closed-form kinetics with blocking (`adhesion_kinetics`) in 1-, 2-
      and 3-D;
    - its Langmuir and jamming limits;
    - each patch against its own column;
    - the ledger;
    - finite differences of the Jacobian.

    Without a substratum, runs are bit-identical to before. See
    docs/validation.md, "Adhesion to surfaces".
  - **Parameters,** graded in docs/parameters.md §7. The one contrast between
    materials that data support, titanium against zirconia, is a calibration.
    Every binding rate is illustrative.

- **pH, and the elements of saliva: Stage S, increment S1, part one**
  (`marse.chemistry`, [docs/networks.md](docs/networks.md#acids-bases-and-ph)).
  A network can now set a pH in every voxel, and rates can depend on it. The
  equations are in docs/theory.md §3.8.
  - **Acids and bases.** A component with `acid_base: {"pka": [...]}` is an
    acid-base total, such as lactic acid and lactate together, written in its
    most protonated form. Protons are never components. The hydrogen ion
    concentration in each voxel is the unique root of the charge balance over
    the totals and the ions of fixed charge, found by a bracketed Newton
    method in log h. `pkw` sets water's ion product.
  - **pH in rates.** A `ph` factor is the cardinal pH model of Rosso et al.
    (1995) at the local pH. Its Jacobian runs through every charged component,
    with dh/dc from the implicit function theorem.
  - **Four more elements.** Formulas may hold P, K, Cl and Na. Each is
    balanced in every network that contains it, and the degree of reduction
    counts them at their valence in phosphate and the salt ions, so those hold
    no electrons.
  - **An absolute tolerance per component.** `absolute_tolerance_mol_per_m3`
    may be an object by component; the rest keep the default.
  - **Outputs.** A run whose network sets a pH records it: a `ph` column in
    `trajectory.csv`, `ph.csv` in space (the substratum and the box), a pH
    field in every ParaView frame, and a summary in the manifest. `marse
    check` prints the pKa values and the pH each starting composition implies.
  - **Verified:**
    - the hydrogen ion concentration against bisection of the charge balance
      (2.4e-12 relative);
    - every voxel left neutral to 1e-12 of the charges present;
    - dh/dc and the rate Jacobian against finite differences;
    - a weak acid and its salt giving back the pH they were made at.

    A network without acids, bases or the new elements runs exactly as
    before: the examples' final digests are unchanged. See docs/validation.md,
    "pH from electroneutrality".

- **The mouth over a site of plaque: Stage S, increment S1, part two**
  (`marse.oral`, `marse.core.reservoir`, [docs/networks.md](docs/networks.md#a-salivary-film-and-the-mouth)).
  A domain in space can now stand under a salivary film, renewed from the
  mouth, instead of under a fixed bulk liquid. The equations are in
  docs/theory.md §4.8 and §9.9.
  - **The film** is the top of the box, closed to the air. Saliva replaces
    each of its voxels at u(z) / l: the film's speed at that height over the
    length of plaque it has crossed (Dawes 1989), with the free-surface
    profile whose shear increment E1 uses.
  - **The mouth** follows Dawes's (1983) model of sugar clearance. Its volume
    grows from the resting volume at a salivary flow that tasting sugar
    raises, and a swallow takes it back without changing its concentrations.
    Secreted saliva moves from resting towards stimulated saliva as the flow
    rises.
  - **Solved together.** The mouth's composition is part of the implicit
    system: its unknowns border the box's, and each linear system is solved
    by a Schur complement on the pool. The S1 prototype showed why. A mouth
    solved apart and corrected after each span conserved to rounding, but its
    answer changed by 0.03 pH with the length of the span.
  - **Two ledgers.** The box's books what crossed into the film. The second
    checks the box and the mouth together, against what was secreted and
    swallowed.
  - **Outputs.** `mouth.csv` gives the mouth's volume, flow, swallows,
    composition and pH at every recorded time. The manifest adds both balances
    and the model's version. `marse check` describes the film and the mouth.
  - **Verified:**
    - the mouth's clearance against Dawes's closed form (3.7e-11 over ten
      swallows);
    - the box and the mouth conserving everything over an hour (5.6e-16);
    - the same pH, within 6.1e-5, whether the mouth runs a minute or 5 s
      ahead of the box;
    - the exchange of a film with its pool against its closed form;
    - the bordered linear system against the Jacobian;
    - replay.

    See docs/validation.md, "The mouth and its film".

- **What is eaten and drunk: Stage S, increment S1, part three**
  (`marse.oral.diet`, [docs/networks.md](docs/networks.md#the-diet)). A
  mouth can now take a diet: intakes listed in order, one at a time, each
  from a start for a duration. The equations are in docs/theory.md §4.9.
  - **A rinse** is taken in at once, held without swallowing, and expelled
    down to the resting volume at its end. A Stephan curve is the plaque's
    response to one.
  - **A drink** is sipped steadily and swallowed as the mouth fills. **A
    food** releases what it holds into the saliva steadily, without liquid,
    as a sweet sucked slowly does. The mouth tastes what they bring, and its
    flow rises.
  - **Mixing.** While an intake lasts, the film is mixed with the mouth's
    liquid (Dibdin 1990), once a second unless the intake says otherwise.
  - **Food left on the teeth.** An intake may leave an amount of a
    particulate component in the film when it ends, over a region of the
    substratum. A process of the network releases what dissolves from it, as
    starchy particles held on the teeth release sugars (Kashket, Zhang and
    Van Houte 1996).
  - **Booked.** The ledger of the box and the mouth counts what was eaten and
    expelled; the box's own counts the food placed in it. The manifest records
    the intakes taken, what was eaten and expelled, and the diet's model
    version. `marse check` lists the diet.
  - **The mouth expects the plaque to go on giving sugar back** at the rate it
    just did, when it runs ahead of the box. After a rinse, the mouth's sugar
    had changed by 1% with the length of the span; it now changes by 8e-4.
    Runs under a mouth change accordingly.
  - **Verified:**
    - every intake booked as eaten within 5e-15 of what was stated, and the
      box and the mouth conserving everything over an hour with a rinse, a
      sipped drink and a sweet that sticks (1.2e-15);
    - the same pH, within 8.6e-5 over an hour after a rinse, whether the
      mouth runs a minute or 5 s ahead of the box;
    - a rinse held and expelled, a drink setting the mouth's sugar, a sweet
      releasing its sugar, and food placed only over its region;
    - a diet that starts after the run changing nothing, bit for bit;
    - replay.

    See docs/validation.md, "The diet".

- **The Stephan curve: Stage S, increment S1, part four**
  (`examples/environments/oral`, `examples/stephan_curve.py`,
  [docs/environments.md](docs/environments.md#the-oral-scenes)). Three oral
  scenes give 150 µm of plaque under the mouth the same sugar in three ways:
  a rinse of 10% sucrose held for a minute, 100 mL of it sipped over 20
  minutes, and the rinse with food left on the teeth. Saliva's buffers are
  those Bardow et al. (2000) measured, at rest and stimulated.
  - **The rinse gives a Stephan curve**, meeting every criterion set before
    the stage was built: the pH falls 1.9 units, to 4.87 at 16 minutes, and
    is back above 6 at 43 minutes, with 15 mM more lactate in the plaque at
    7 minutes. A two-hour curve runs in 13 s.
  - **What keeps plaque acid for longer.** Sipping keeps the plaque below pH
    5.5 for 52 minutes and food left on the teeth for 56, against 32 after
    the rinse. Low salivary flow, a slower film, thicker plaque and more
    fixed buffer each move the curve the ways the literature reports, and
    each is a test.
  - **Fixed charges release their counter-ions.** A plaque buffer holding its
    cations fixed turned the whole salivary film above it to pH 3.3 after a
    rinse, as the lactate leaving the plaque carried the acid's protons with
    it. Two fast processes now keep the cations the groups hold equal to
    their charge, so the acid stays on the buffer until the saliva's alkali
    takes it off, and the film is never more acid than the plaque. A new
    rate factor, `dissociated`, gives the protons an acid-base total has lost
    at the local pH, with its Jacobian through dh/dc. The scenes give every
    charged component one diffusivity, so that diffusion separates no charge.
  - **Calibrated again.** With the acid held in the plaque, the prototype's
    300 µm of plaque stayed acid for more than an hour. The plaque's
    thickness, the film's speed, the buffer and the ionic diffusivity were
    calibrated to the criteria, within their ranges, all confidence C
    ([parameters.md §8](docs/parameters.md#8-saliva-plaque-and-diet)).
  - **The pH solve can no longer cycle.** In one voxel of these scenes,
    Newton's method stepped for ever between the two ends of its bracket. A
    step must now halve the one before, or the bracket is bisected (Press et
    al. 2007); 2,000 plaque-like voxels take 12 iterations instead of 19.
  - `marse.oral.stephan` measures a curve: its minimum and when, the minutes
    and the area below pH 5.5, and when it is back above 6. `marse check`
    says that a column is solved directly, and calls its fastest time scale
    a process rather than growth.

    See docs/validation.md, "The Stephan curve".

- **Plaque that spreads, wears and is brushed off: Stage S, increment S2,
  part one** (`marse.biofilm.spreading`,
  [docs/networks.md](docs/networks.md#plaque-that-spreads)). A column's
  plaque now has room to grow. A domain may give it a `plaque`:
  - the packing concentration of each component that takes up space;
  - the components carried with them;
  - a maximum height;
  - a wear velocity.

  The equations are in docs/theory.md §6.1 and §6.3.
  - **Spreading.** After every step, the solid is packed up the column by its
    cumulative volume (`displacement_1d_v1`; Wanner and Gujer 1986). The
    column is left with full voxels, then at most one partly filled. Every
    component moves in the proportions of the voxel it came from. The packing
    conserves exactly and makes nothing negative. Spreading in two and three
    dimensions stays with Stage 2d, because there the choice of mechanism
    changes conclusions.
  - **Detachment and wear.** Solid pushed past the maximum height leaves the
    column. A stated wear velocity takes the surface off: half before each
    step and half after, which keeps the step second order. Without a mouth,
    what leaves goes to the bulk liquid. Under the mouth, it enters the
    mouth's saliva and is swallowed.
  - **Brushing and flossing.** A domain's `hygiene` lists timed cleanings.
    Each takes a share of the plaque off from its surface down, and the same
    share of any food left on the teeth. A brushing takes 42% unless told
    otherwise (Slot et al. 2012). A flossing must state its share. The mouth
    expels what a cleaning takes, and both ledgers book it.
  - **Outputs.**
    - `plaque.csv` records the plaque's thickness, and what each filling
      component holds, has detached and has been brushed off.
    - The manifest records the spreading provider and what the plaque lost.
    - `marse check` describes the plaque and its cleanings.
  - **Verified:**
    - against the closed form for a film growing against wear, to 9.7e-7;
    - against the packing done voxel by voxel, on 500 random columns;
    - with a brushing under the mouth that leaves exactly 58% of the plaque,
      with both ledgers closed to 1.9e-16.

    A domain without a plaque runs as before. See docs/validation.md,
    "Plaque that spreads".

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

### Fixed

- **The mouth's books close to rounding on long runs.** Each span under the
  mouth was integrated on the run's clock, so a step of a few seconds was the
  difference of two times hours into the run, which loses digits in
  proportion to the time. What the steps added and what the span booked
  drifted apart: over six hours of meals, the box and the mouth together
  balanced to only 4e-13 of the sugar eaten, and a day or a year would have
  drifted further. Each span now keeps its own clock, from zero, and the
  same six hours balance to 3e-16. **Results changed**, in their last digits:
  the three oral scenes give the same Stephan curves to six figures, with new
  final digests. Found while building Stage S2's day-long scenes.
- **Runs in space now replay bit for bit on any number of threads.** They used
  to differ in their last digits between machines with different numbers of
  cores. The coarsest multigrid level was inverted by LAPACK, and OpenBLAS
  factorises a matrix of more than about 100 unknowns in parallel, so the
  result depended on its thread count. The 1-D example's final digest was
  different at 1, 2 and 4 threads, and `marse replay` reported a difference
  on a machine with a different core count.
  - **The fix.** MARSE now inverts that level itself, by a blocked LU
    factorisation with partial pivoting, which agrees with LAPACK to about one
    unit in the last place. Its matrix products run through OpenBLAS, which
    computes them the same way on any number of threads. The rest is numpy's
    elementwise arithmetic.
  - **Faster assembly.** The level's matrix is assembled with one operator
    application per voxel instead of one per unknown. It is the same matrix,
    bit for bit.
  - **Tested.** A test runs the 1-D example in subprocesses at 1 and 4
    threads and requires one final digest.
  - **Results changed.** Every run in space differs from before in its last
    digits, and no more.
- `marse` no longer switches every NumPy floating-point error to "warn" when a
  command ends. It now restores the caller's settings, where before it turned
  underflow into a warning, or into an error under `pytest`, for the rest of the
  process.
