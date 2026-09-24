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
