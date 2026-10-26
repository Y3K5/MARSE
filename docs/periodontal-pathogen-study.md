# Periodontal pathogen biofilm study

**Not a calibrated pocket.** The study to run is the version-2 column
`examples/networks/periodontal_pathogen_v2.json`. Three obligate anaerobes
share alanine, a stand-in for a peptide pool, in one basal layer of a film.
Carbon, nitrogen and electrons balance. Oxygen is not a growth substrate.
The yields, rates and the generic biomass formula are placeholders. The
column is not the old 60 by 40 cell map, and two-dimensional biomass
spreading is not implemented, so the populations do not crawl sideways.
Nothing here is a pocket, a virulence mechanism, a patient outcome, or a
strain-specific physiology.

The version-1 file `periodontal_pathogen_biofilm.json` is still runnable with
`marse ecosystem`. That engine is unverified. Its defects are listed in
[validation.md](validation.md#the-two-dimensional-ecosystem-engine-is-not-yet-verified),
and the command prints a warning. Those numbers are not a result to cite.
The version-2 column does not inherit that engine's mass-balance, oxygen,
mutation or carrying-capacity behaviour. It also does not retire the
version-1 defects inside `marse ecosystem`; that code is unchanged.

`marse bench` scores the column's carbon, nitrogen and electron balance,
and records the version-1 file only as an unverified engine. The two files
are not compared as biomass. The cases, the limit of that score, and the
plan for retiring version 1 without deleting it are in
[periodontal.md](benchmarks/periodontal.md).

## What the version-2 column connects

- **One shared substrate.** Alanine limits all three populations. There is no
  competition matrix and no per-species carrying capacity.
- **Shared room in the column.** Continuum spreading moves biomass that no
  longer fits, and the three populations count toward the same voxel.
- **Anaerobic growth.** Each organism is an obligate anaerobe. No process
  consumes or requires oxygen.
- **Not carried over from version 1:** mutation of grid cells, cardinal
  temperature, a moisture scalar, adhesion edges, and lateral crawling.

## What the version-1 file connected

This is the structure of `periodontal_pathogen_biofilm.json` only.

- **Shared-resource competition:** all three species use oxygen, heme,
  peptides, and carbon, with an explicit competition matrix.
- **Spatial overlap:** each species starts in a separate basal patch and
  spreads through the same grid.
- **Attachment and detachment:** biomass is biased toward the bottom surface,
  with species-specific adhesion and detachment rates.
- **Declared production fields:** growth contributes small, named nutrient
  production terms as bookkeeping proxies for cross-feeding or released
  resources. These are not assigned to a particular metabolite.
- **Named metabolite proxies:** `succinate` and `acetate` are transported and
  decayed as explicit additive fields with declared producer/consumer rates.
  Their names do not imply calibrated exchange stoichiometry.
- **Mutation flags:** low-probability growth multipliers provide a way to
  visualize sensitivity to heritable variation, not a gene-level model.

The first analyses should compare single-species controls against the
three-species run: total biomass, occupied area, overlap, resource depletion,
and which species persists in shared regions. A useful next extension is to
replace the production proxies with evidence-linked metabolites and to add
experimentally measured growth and competition coefficients.

## Run

```bash
marse check examples/networks/periodontal_pathogen_v2.json
marse run examples/networks/periodontal_pathogen_v2.json -o runs/periodontal-v2
marse replay runs/periodontal-v2/manifest.json
```

The run writes the usual spatial files plus `summary.json` and `report.html`.
Open `runs/periodontal-v2/vtk/run.pvd` in ParaView for the column.

The version-1 viewer is still available, and still unverified:

```bash
marse ecosystem examples/experiments/periodontal_pathogen_biofilm.json \
  -o runs/periodontal-pathogen-biofilm
```

Use **Fields + particles** there to see continuous biomass fields with
deterministic population representatives. The particles are visualization
aids. That simulation remains the version-1 spatial biomass-field model.

## Controls and environmental variants

Run the compact control and sensitivity workflow without exporting every frame:

```bash
python examples/periodontal_variant_analysis.py \
  --output runs/periodontal-variant-analysis
```

It writes `summary.csv` and `manifest.json`, not large trajectory files. The
workflow includes the mixed biofilm and one single-species control for each
organism at 30, 37, and 40 C, with three deterministic replicate seeds at
reference and reduced moisture. The controls help separate shared-resource
competition from species-specific growth.

Temperature is applied through each species' declared cardinal capability.
Moisture is currently an explicit dimensionless condition whose illustrative
reduced level scales available oxygen and carbon to 65% of reference. This is
a sensitivity proxy, not a validated periodontal moisture law; replace it
with measured medium or pocket data before quantitative conclusions.

The manifest records the base configuration checksum, software environment,
temperature levels, moisture levels, replicate set, and scenario count. Keep
generated runs outside Git and archive a selected release plus the manifest
with a DOI when preparing a paper.

## Parameter provenance and calibration status

`periodontal_pathogen_provenance.json` is a sidecar register for the
experiment. It records parameter paths, units, the current source identifier,
evidence grade, calibration status, and the measurements needed to replace
placeholders. Its current status is explicitly `exploratory-placeholder`:
the repository does not yet claim organism-specific measurements for these
values.

The variant workflow hashes both the experiment and provenance sidecar into
`manifest.json`. Before quantitative publication, replace the placeholder
source with source-linked growth curves, uptake/yield data, pairwise and
three-species competition observations, attachment/detachment measurements,
and named metabolite measurements. Calibration should then be performed
against controls before interpreting mixed-community outcomes.
