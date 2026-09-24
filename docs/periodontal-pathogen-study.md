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
