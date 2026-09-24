# Scientific culture evidence

`marse.science` stores the conditions and methods attached to culture
observations. It is deliberately not a table of universal optima: a result is
only transferable when its organism or strain, medium, physical format,
atmosphere, inoculum, incubation, endpoint, and source context are compatible.

## Evidence records

- `SourceRecord` identifies a paper, standard, or curated strain database and
  carries the conservative A/B/C evidence grade used in
  [`parameters.md`](parameters.md).
- `MediumRecipe` names a versioned composition and links it to sources.
- `AgarProtocol` records gel concentration, plate depth, volume, inoculation,
  and atmosphere. Agar concentration is a physical transport parameter, not a
  species-independent growth optimum.
- `MeasurementMethod` identifies the endpoint and unit. Optical-density
  measurements require a condition-specific calibration record.
- `KineticObservation` stores fitted growth results separately from the raw
  culture context, including model family, fitting interval, replicates, and
  uncertainty.
- `CultureRecord` joins those records for one strain-specific observation.

Datasets are loaded from JSON with `load_culture_dataset` and serialized in
canonical form for a reproducible SHA-256 checksum. Missing references,
unsupported endpoints, invalid ranges, and universal optimum claims are
rejected before a record can be used to construct a simulation.

## Evidence-to-experiment compilation

`compile_culture` is the first bridge from evidence to execution. It compiles
a strain-specific **broth** record with a fitted growth rate into a validated
`ExperimentConfig`, while requiring the experiment designer to provide
substrate concentration, yield, half-saturation, initial biomass, duration,
and timestep. Those values are not safely inferable from a culture record.

The returned `ExperimentCompilation` includes the selected sources, dataset
checksum, assumptions, and unresolved requirements. Agar, semi-solid, and
biofilm records are intentionally not projected into the well-mixed batch
kernel; they remain evidence until a spatial compiler is implemented.

## Growth-curve calibration

`marse.calibration` accepts strictly increasing time points and one or more
positive replicate trajectories. It supports exponential, logistic, Gompertz,
and Baranyi screening fits and returns parameters, pointwise residuals, RMSE,
replicate count, and identifiability warnings. Log transforms are explicit and
must match the endpoint's measurement scale.

These fits are deterministic calibration aids. They do not estimate confidence
intervals, correct for censoring, or establish causal biological parameters.
Those require an experiment-specific statistical model and should be added only
when the observation design supplies enough information.
