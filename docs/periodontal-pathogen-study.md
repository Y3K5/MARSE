# Periodontal pathogen biofilm study

The `periodontal_pathogen_biofilm.json` experiment is a reproducible,
exploratory starting point for comparing *Treponema denticola* (T. denticola),
*Tannerella forsythia* (T. forsythia), and *Porphyromonas gingivalis*
(P. gingivalis) in one spatial biofilm-like environment.

It does **not** claim to reproduce a periodontal pocket, virulence mechanism,
patient outcome, or strain-specific physiology. The numeric values are
illustrative placeholders that make the interaction structure executable.
They should be replaced with organism- and medium-specific measurements before
quantitative interpretation.

## What the model connects

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

## Run and visualize

```bash
python -m marse.cli ecosystem \
  examples/experiments/periodontal_pathogen_biofilm.json \
  --output runs/periodontal-pathogen-biofilm
open runs/periodontal-pathogen-biofilm/viewer.html
```

Use **Fields + particles** to see continuous biomass fields with deterministic
population representatives. The particles are visualization aids; the
underlying simulation remains a spatial biomass-field model.

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
