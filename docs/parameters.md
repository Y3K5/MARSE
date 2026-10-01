# Parameter reference

Numerical inputs for the equations in [`docs/theory.md`](theory.md), with
units, sources and an explicit confidence level for each.

This is a **modelling parameter reference**, not a laboratory protocol. It
records where published rate constants came from and under what conditions
they were measured, so that a simulation can state the provenance of every
number it used.

## How to read this table

MARSE's design rule is that every parameter carries units, a source and a
confidence note (see [`docs/specification.md`](specification.md)). The
verification column here is deliberately conservative:

| Mark | Meaning |
|---|---|
| **A** | Value computed from a published closed-form correlation implemented in MARSE and checked against independent reference values in the test suite. |
| **B** | Bibliographic record confirmed (authors, journal, volume, pages, DOI); the numerical value is consistent with multiple secondary sources but the primary text was **not** opened during preparation. |
| **C** | Illustrative order-of-magnitude default. Must be replaced with a fitted or cited value before use in any published analysis. |

Only **A** values should be treated as settled. **B** values need checking
against the primary source before they appear in a paper; they are recorded
here so the check has a defined target. **C** values exist to make examples
runnable and are labelled as such in configuration.

> The network available while preparing this document blocked direct access to
> publisher sites, PubMed and Crossref, so no value below is marked **A**
> unless MARSE computes it from a formula it implements and tests. This is a
> statement about the preparation of the document, not about the quality of
> the underlying literature.

---

## 1. Physical constants and correlations

These are computed by MARSE from published correlations and verified in
`tests/test_kinetics.py`.

| Quantity | Value | Unit | Source | Verified |
|---|---|---|---|---|
| O₂ solubility, fresh water, air, 1 atm, 20 °C | 9.09 | mg L⁻¹ | Benson & Krause (1984) | **A** |
| O₂ solubility, 25 °C | 8.26 | mg L⁻¹ | Benson & Krause (1984) | **A** |
| O₂ solubility, 30 °C | 7.56 | mg L⁻¹ | Benson & Krause (1984) | **A** |
| O₂ solubility, 37 °C | 6.73 (= 0.210 mM) | mg L⁻¹ | Benson & Krause (1984) | **A** |
| O₂ diffusivity in water, 25 °C | 2.00 × 10⁻⁹ (2000 µm² s⁻¹) | m² s⁻¹ | Han & Bartels (1996) | **A** |
| O₂ diffusivity in water, 37 °C | 2.62 × 10⁻⁹ (2624 µm² s⁻¹) | m² s⁻¹ | Han & Bartels (1996) | **A** |
| O₂ molar mass | 31.998 | g mol⁻¹ | definition | **A** |
| Water viscosity, 25 °C / 37 °C | 0.890 / 0.692 | mPa s | Huber et al. (2009), IAPWS 2008 | **B** |
| Glucose diffusivity in water, 25 °C | ≈ 600 | µm² s⁻¹ | compiled value; see note | **B** |
| Glucose diffusivity, 37 °C (Stokes–Einstein scaled) | ≈ 800 | µm² s⁻¹ | derived, theory.md §4.6 | **C** |
| CO₂ diffusivity in water, 25 °C | ≈ 1.9 × 10⁻⁹ | m² s⁻¹ | commonly tabulated; primary source not yet checked | **C** |
| NH₄⁺ diffusivity in water, 25 °C | ≈ 2.0 × 10⁻⁹ | m² s⁻¹ | commonly tabulated limiting ionic value; primary source not yet checked | **C** |
| Lactate diffusivity in water, 25 °C | ≈ 1.0 × 10⁻⁹ | m² s⁻¹ | order of magnitude; primary source not yet checked | **C** |

**Diffusivities in the spatial example.** `examples/networks/surface_biofilm_3d.json`
starts from these values. It scales them to 37 °C by Stokes–Einstein (theory.md
§4.6; a factor of about 1.34 from 25 °C) and multiplies them by the class means
of §2 below, giving effective values inside a biofilm, in m² s⁻¹:

| Solute | Value |
|---|---|
| oxygen | 1.13 × 10⁻⁹ |
| glucose | 2.3 × 10⁻¹⁰ |
| ammonium | 1.47 × 10⁻⁹ |
| carbon dioxide | 1.10 × 10⁻⁹ |
| lactate | 3.9 × 10⁻¹⁰ |

Configuration schema version 2 gives each component one diffusivity for the
whole box. Until diffusivity depends on the biofilm (Stage 3), the liquid above
the colonies diffuses at these effective values too, which understates
transport there by a factor of 2 to 3. The CO₂, ammonium and lactate values
are graded C and must be checked against a primary source before any published
analysis.

Implemented in `marse.spatial.solutes`. Both oxygen correlations refuse to
extrapolate outside their stated validity ranges (0–40 °C and 0–95 °C
respectively) rather than returning a plausible-looking wrong number.

**Corrections MARSE does not apply automatically.** Dissolved salts lower
oxygen solubility, so a culture medium holds less than pure water; a 5 % CO₂
atmosphere displaces oxygen to roughly 0.95 × the tabulated value. Both must
be stated in configuration.

### Composition used by reaction networks

Reaction networks ([`docs/networks.md`](networks.md)) balance carbon, nitrogen
and electrons from each component's formula alone. The atomic weights only
convert amounts to masses, for the molar masses `marse check` prints.

| Quantity | Value | Unit | Source | Verified |
|---|---|---|---|---|
| Atomic weights of C, H, N, O | 12.011, 1.008, 14.007, 15.999 | g mol⁻¹ | definition: IUPAC abridged standard atomic weights (Prohaska et al. 2022) | **A** |
| Average elemental composition of microbial biomass | CH₁.₈O₀.₅N₀.₂: 4.2 e⁻ and 24.63 g (ash-free) per C-mol | per C-mol | Heijnen & van Dijken (1992); Roels (1983) | **B** |

The biomass composition is an average. Real biomass varies between organisms
and with growth conditions, so a network states each organism's formula
explicitly and MARSE supplies no default. The molar masses are checked in
`tests/test_schema_formula.py`, glucose at 180.156 g mol⁻¹.

---

## 2. Effective diffusivity in biofilms

Relative effective diffusivity $f = D_e/D_{aq}$ (theory.md §4.3), grouped by
solute physical chemistry rather than by organism.

| Solute class | Mean $f$ | Source | Verified |
|---|---|---|---|
| Inorganic ions | ≈ 0.56 | Stewart (1998) | **B** |
| Small nonpolar solutes, MW ≤ 44 (includes O₂) | ≈ 0.43 | Stewart (1998) | **B** |
| Organic solutes, MW > 44 (includes glucose) | ≈ 0.29 | Stewart (1998) | **B** |

These are **means of a wide distribution**, not constants. $f$ varies within a
single biofilm, changes as the biofilm matures, and depends on density and
matrix composition. Any conclusion sensitive to gradients should report a
sensitivity analysis over $f$, not a single value.

---

## 3. Growth rates of model organisms

Reference values for *Escherichia coli* K-12, the best-characterised
laboratory model organism and MARSE's primary calibration target.

| Quantity | Value | Unit | Conditions | Source | Verified |
|---|---|---|---|---|---|
| Doubling time, rich medium, 37 °C | ≈ 20–30 | min | LB, aerated | Sezonov et al. (2007) | **B** |
| → equivalent $\mu$ | 1.4–2.1 | h⁻¹ | derived, $\ln 2/t_d$ | — | **A** |
| Doubling time, minimal medium + glucose, 37 °C | ≈ 38–67 | min | M9/N⁻C⁻, strain-dependent | compiled; Reshes et al. (2008), Soupene et al. (2003) | **B** |
| → equivalent $\mu$ | 0.62–1.09 | h⁻¹ | derived | — | **A** |
| Biomass yield on glucose, aerobic | ≈ 0.45 | g dry wt (g glucose)⁻¹ | batch, minimal medium | compiled (BioNumbers); Link et al. (2008) | **B** |
| Specific O₂ uptake rate | ≈ 20–27 | mmol O₂ (g dry wt)⁻¹ h⁻¹ | aerobic, glucose | compiled (BioNumbers) | **B** |
| $K_S$, glucose | µg L⁻¹ to mg L⁻¹ range | — | strongly method-dependent | Kovárová-Kovar & Egli (1998) | **B** |

**On $K_S$ specifically.** Kovárová-Kovar & Egli (1998) document that reported
half-saturation constants for the same organism and substrate span **orders of
magnitude**, depending on cultivation history and measurement method. $K_S$ is
not a fixed property of a species. Treat any single literature value as one
sample from a broad distribution, and propagate that uncertainty rather than
reporting a point estimate. This is the largest single source of parameter
uncertainty in the models of theory.md §1.3.

### Consistency check

A useful arithmetic check that these numbers hang together: 0.4 % (w/v)
glucose is 22.2 mM, and at a yield of 0.45 g g⁻¹ it supports

$$
\Delta X = 0.45 \times 4\ \text{g L}^{-1} = 1.8\ \text{g L}^{-1}
$$

of dry biomass — a realistic final density for a minimal-medium batch culture.
A parameter set that fails this kind of check is wrong regardless of its
provenance.

---

## 4. Cardinal temperatures and pH

For the secondary models of theory.md §2. Cardinal values are
**strain- and condition-dependent**, and published estimates for the same
species differ by several degrees.

### 4.1 Temperature (CTMI)

| Organism | $T_{\min}$ | $T_{\mathrm{opt}}$ | $T_{\max}$ | Source | Verified |
|---|---|---|---|---|---|
| *E. coli* K-12 | ≈ 5–8.5 °C | ≈ 40 °C | ≈ 46–47 °C | Van Derlinden et al. (2008); Baka et al. (2013) | **B** |

Reported behaviour that constrains $T_{\max}$: growth curves are smooth at
40–43 °C, the exponential phase is interrupted at 44–45 °C, and at 46 °C a
period of growth is followed by inactivation. The CTMI requires
$T_{\mathrm{opt}} \ge \tfrac{1}{2}(T_{\min} + T_{\max})$ (theory.md §2.2), which
these values satisfy.

MARSE's examples use $(T_{\min}, T_{\mathrm{opt}}, T_{\max}) = (6, 40, 47)$ °C
as an **illustrative** set consistent with that range, marked **C** in
configuration. Fitting cardinal values to the specific strain and medium under
study is part of Phase 3, not a lookup.

### 4.2 pH (CPM)

| Organism | $\mathrm{pH}_{\min}$ | $\mathrm{pH}_{\mathrm{opt}}$ | $\mathrm{pH}_{\max}$ | Source | Verified |
|---|---|---|---|---|---|
| *E. coli* | ≈ 4 | ≈ 7 | ≈ 9–10 | Presser et al. (1997); Rosso et al. (1995) | **B** |

Presser et al. (1997) report growth at pH 4.0 but not at pH 3.7 in the absence
of organic acid, which brackets $\mathrm{pH}_{\min}$ between those values.

**Undissociated organic acids inhibit far more strongly than pH alone**, so
a pH-only $\gamma$ factor will over-predict growth in an acidified,
fermenting environment. Where organic acids matter, they belong in the model
as diffusing solutes with their own inhibition term, not folded into
$\mathrm{pH}_{\min}$.

---

## 5. Biofilm structure

| Quantity | Value | Conditions | Source | Verified |
|---|---|---|---|---|
| O₂ penetration depth | tens of µm | dense colony biofilm | Walters et al. (2003) | **B** |
| Zone of active protein synthesis | ≈ 30–60 µm from the oxic surface | colony and flow-cell biofilms | Werner et al. (2004) | **B** |
| Total biofilm thickness | often several × the penetrated depth | — | as above | **B** |

The relationship between these three rows *is* the stratification result of
theory.md §5.1: an active surface layer of tens of micrometres over a much
thicker inactive interior. MARSE's validation case V3 targets this
qualitative structure — the existence and approximate scale of the active
zone — rather than a specific micrometre value, because that is what the
underlying physics predicts robustly.

### Derived consistency check

Using $D_e = 0.43 \times 2624 = 1128$ µm² s⁻¹ and $C_0 = 0.210$ mM at 37 °C,
a 50 µm penetration depth requires a volumetric uptake of 0.19 mM s⁻¹
(theory.md §5.2), which at 20–27 mmol gDW⁻¹ h⁻¹ implies roughly
25 g L⁻¹ dry biomass. Independent parameters from three different sources
produce a mutually consistent picture, which is weak evidence that none is
badly wrong.

---

## 6. Adaptation and phenotypic switching

| Quantity | Value | Source | Verified |
|---|---|---|---|
| Phenotypic switching rates (persistence) | measurable, low, strain-dependent | Balaban et al. (2004) | **B** |

Balaban et al. (2004) established that switching between normally growing and
slow-growing states occurs at quantifiable rates and distinguished
phenotypes that arise on transition to stationary phase from those generated
continuously during growth. MARSE's transition-rate parameters (theory.md §8)
take this form: rates with units, sources and confidence, not tuned constants.

No default switching rates are supplied. A transition rate with no source is
a fitted parameter and must be declared as one.

---

## 7. Surfaces and adhesion

Cells binding to a surface (theory.md §6.4) need two kinds of number:

- **Delivery.** MARSE computes how fast cells reach the surface from physics.
- **Binding.** How many bind, and how firmly, it cannot compute. Those numbers
  are stated in each scene.

For the laboratory and dental scenes of [`docs/environments.md`](environments.md),
most binding parameters are illustrative. The literature supports one contrast
between materials, titanium against zirconia. This section records which
numbers are which.

### 7.1 Physics and measurements

| Quantity | Value | Unit | Source | Verified |
|---|---|---|---|---|
| Lévêque constant, $1/(\Gamma(4/3)\,9^{1/3})$ | 0.5384 | — | Lévêque (1928); matches a numerical solution of the boundary-layer problem to 1% | **A** |
| Brownian diffusivity of a 1 µm sphere in water, 25 °C | 0.4907 | µm² s⁻¹ | Stokes–Einstein, with the water viscosity of §1 | **A** |
| Jamming limit of random sequential adsorption of disks | 0.547 | area fraction | Feder (1980) | **B** |
| Initial deposition rates of oral streptococci onto glass in a parallel-plate flow cell | 0 to 2.9 × 10³ | cm⁻² s⁻¹ | Sjollema, Busscher & Weerkamp (1988) | **B** |
| Adhesion force of streptococci on saliva-coated enamel, contact of 0 → 120 s | −0.7 → −10.3 | nN | Mei et al. (2009) | **B** |
| Thickness of the salivary film | 70–100 | µm | Collins & Dawes (1987) | **B** |
| Velocity of the salivary film, by site in the mouth | 0.8–7.6 | mm min⁻¹ | Dawes et al. (1989) | **B** |
| Wall shear rate of the salivary film, $3\bar u/\delta$ | 0.4–5.4 | s⁻¹ | derived from the two rows above | **A** (formula) |
| Threshold roughness for plaque retention | $R_a$ ≈ 0.2 | µm | Bollen, Lambrechts & Quirynen (1997) | **B** |

The roughness threshold is recorded for when retention is modelled. Binding
does not use it, because roughness shelters cells from removal more than it
changes how they bind.

### 7.2 Titanium against zirconia

Four studies bear on the one contrast between materials that the dental scene
draws:

- **In the mouth.** After 24 h on discs worn by volunteers, bacteria covered
  19.3 ± 2.9% of titanium and 12.1 ± 2.0% of zirconium oxide (Scarano et al.
  2004; **B**). The ratio is 0.63.
- **Early colonizers in the laboratory.** *S. sanguinis*, *S. gordonii* and
  *S. oralis* adhered significantly less to polished zirconia than to polished
  titanium, while *S. mutans* adhered equally to both (Oda et al. 2020; **B**).
- **Surface energy.** Polished zirconia had lower surface free energy than
  polished titanium, and less *S. mitis* and *Prevotella nigrescens* adhered
  to it, especially under a salivary pellicle. The authors name surface free
  energy as the main factor on smooth surfaces (Al-Radha et al. 2012; **B**).
- **A study that found no difference.** On low-roughness yttria-stabilized
  zirconia and titanium carried in the mouth, initial adhesion and biofilm
  formation were comparable (Al-Ahmad et al. 2016; **B**).

The dental scene therefore sets zirconia's attachment efficiency to 0.63 of
titanium's, calibrated on Scarano et al. The 24-hour ratio the scene produces
reproduces that calibration. It is not a prediction, and Al-Ahmad et al.
(2016) is a reason to doubt it.

No comparison was found for enamel or acrylic (PMMA) against titanium. The
scene gives all three the same efficiency, so any difference between them in
its output would be a bug.

### 7.3 Illustrative values in the scenes

Every value in this table is **C**. Each one makes the scenes run with
plausible magnitudes and must be replaced by a measured value before use in
any published analysis.

| Quantity | Lab scene | Dental scene | Unit | Basis |
|---|---|---|---|---|
| Attachment efficiency $\alpha$ | 0.5 bare glass, 0.3 saliva-coated | 0.5 enamel, titanium and PMMA; 0.315 zirconia | — | Keeps the lab scene's initial deposition inside the range of §7.1. Zirconia is 0.63 × titanium (§7.2). |
| Detachment $k_{\mathrm{off}}$ | 30 bare, 10 coated | 10 | h⁻¹ | Reversibly bound cells leave within minutes. |
| Locking $k_{\mathrm{lock}}$ | 90 | 90 | h⁻¹ | A time constant of 40 s, matching the two-minute strengthening of §7.1. |
| Cells in the bulk liquid | 3 × 10⁸ *S. oralis* | 1 × 10⁷ *S. oralis*, 5 × 10⁶ *S. sanguinis* | mL⁻¹ | A dense flow-chamber suspension, and a share of saliva's order of 10⁸ bacteria per mL. |
| Cell diameter | 0.9 | 0.9 *S. oralis*, 1.0 *S. sanguinis* | µm | Streptococci are cocci about a micrometre across. |
| Carbon per cell | 7 | 7 *S. oralis*, 9 *S. sanguinis* | fmol | A cell of this size holds a few to ten femtomoles of carbon. |
| Blocked area per bound cell | 10 | 10 | µm² | 13 to 16 times a cell's footprint. It gives a jamming density of 5.5 × 10⁶ cm⁻². |
| Liquid | water at 37 °C, 0.692 mPa s (§1) | saliva at 35 °C, 1.2 mPa s | — | Whole saliva is viscoelastic and shear-thinning, and its viscosity depends on stimulation. The transfer velocity scales as $\eta^{-2/3}$. |
| Wall shear rate | 15 | 1.25: a film 80 µm thick moving at 2 mm min⁻¹ | s⁻¹ | A typical flow-chamber shear, and a film inside the ranges of §7.1. |
| Distance downstream | 20 | 5 | mm | Where on the chamber's plate, or the tooth, the scene sits. |
| Growth | none (buffer) | μmax 0.3 and 0.25 h⁻¹; $K_S$ 0.1 mM glucose, 0.01 mM ammonium; yield 0.6 C-mol per mol glucose | — | Lactic fermentation of glucose, as an order of magnitude. |
| Salivary glucose, ammonium, CO₂ | — | 0.2, 2 and 1 | mM | A daytime average; meals arrive with the saliva-and-diet stage. |

The transfer velocity depends only weakly on what is uncertain. It scales as
the diffusivity to the power 2/3, so as $d^{-2/3}$ and $\eta^{-2/3}$, and as
the shear to the power 1/3. If saliva at these low shear rates were ten times
as viscous as the value used, cells would arrive 4.6 times more slowly. The
number bound scales directly with $\alpha$ and the suspension, which are the
least certain of all.

---

## 8. Saliva, plaque and diet

The oral scenes of [`docs/environments.md`](environments.md#the-oral-scenes) run
a column of plaque under a salivary film, renewed from a mouth that secretes,
swallows and takes a diet (theory.md §3.8, §4.8 and §4.9). The film's surface
may be at the air (theory.md §4.10), which holds oxygen at its solubility in
§1 above. Saliva's buffers and the mouth's volumes are measured values; what
the plaque does with sugar is calibrated, all of it **C**, to the criteria for
a Stephan curve set before Stage S1 was built (docs/validation.md, "The
Stephan curve").

### 8.1 Measured values

| Quantity | Value | Unit | Source | Verified |
|---|---|---|---|---|
| Resting whole saliva: pH; bicarbonate; phosphate; flow | 6.8; 4.4 (5.3 with its CO₂, at 29.3 mmHg); 4.5; 0.55 | —; mM; mM; mL min⁻¹ | Bardow et al. (2000), collected under oil so that no CO₂ was lost | **B** |
| Stimulated whole saliva: pH; bicarbonate; phosphate; flow | 7.2; 9.7 (10.5 with its CO₂, at 25.7 mmHg); 3.8; 1.66 | —; mM; mM; mL min⁻¹ | Bardow et al. (2000) | **B** |
| Conditional pKa of carbonic acid; second pKa of phosphate | 6.1; 7.2 | — | Implied by the same saliva: its pH, bicarbonate and CO₂ give 6.11 at rest and stimulated, and its HPO₄²⁻ fractions give 7.19 and 7.20 | **B** (derived) |
| Volume of saliva in the mouth after a swallow, and before one | 0.77 and 1.07 | mL | Lagerlöf & Dawes (1984) | **B** |
| Unstimulated salivary flow | 0.3 | mL min⁻¹ | Dawes (1983) models clearance at about this flow; Bardow et al. (2000) measured 0.55 | **B** |
| Thickness of the salivary film; its velocity by site | 70–100; 0.8–7.6 | µm; mm min⁻¹ | Collins & Dawes (1987); Dawes et al. (1989) | **B** |
| Plaque's buffering | strong at pH 4–5.5, weak near neutrality; about 90% from cell walls and matrix | — | Shellis & Dibdin (1988) | **B** (qualitative) |
| Food retained on the teeth | starchy particles for up to 20 min, accumulating sugars and acids | — | Kashket, Zhang & Van Houte (1996) | **B** (qualitative) |
| Plaque removed by one brushing with a manual toothbrush: overall; by plaque index | 42; 30 to 53 | % | Slot et al. (2012), 59 papers and 212 brushing exercises. MARSE's default for a brushing (theory.md §6.3). | **B** |
| Salivary flow while chewing gum | 10 to 12 times the unstimulated 0.47 in the first minute with flavoured gum; within about 10 min, the flow with gum base alone; still above the unstimulated flow after 2 h | mL min⁻¹ | Dawes & Macpherson (1992); Dawes & Kubieniec (2004). A mouth's `chewing_flow_ml_per_min` is the lasting part (theory.md §4.9). | **B** |

### 8.2 Values in the oral scenes

Every value in this table is **C**. The calibrated ones were chosen so that a
rinse of 10% sucrose meets all five criteria for a Stephan curve with a margin
on each; the others are orders of magnitude.

| Quantity | Value | Unit | Basis |
|---|---|---|---|
| Acid production | 0.5 per C-mol of an 800 C-mM population, so 400 mM of hexose per hour at pH 7 | h⁻¹ | Calibrated. Lactate in the plaque rises by 15 mM at 7 minutes. |
| Half-saturation of acid production for sugar | 1 | mM | An order of magnitude: the plaque ferments at its full rate at the sugar a rinse leaves. |
| Cardinal pH of acid production | 4, 7, 9 | — | One population for all; Stage S2 separates growth from acid production, species by species. |
| Fixed buffer: carboxyl groups; their pKa | 160; 4.8 | mM; — | Calibrated, with its pKa inside the region of strong buffering of §8.1. |
| Exchange of the fixed groups' potassium | 3600 | h⁻¹ | Fast against everything else, so that the bound potassium equals the groups' charge. |
| Thickness of plaque | 150 | µm | Calibrated: sugar lingers in plaque for about three times $L^2/D$, so the thickness sets when the minimum falls. |
| Velocity of the film; the plaque it has crossed | 6; 6 | mm min⁻¹; mm | Calibrated, inside the range of §8.1; the film is renewed every minute on average. |
| Area of plaque the column stands for | 2 | cm² | An order of magnitude for plaque-covered surfaces. |
| Diffusivity of sugar; of every charged component | 2.3 × 10⁻¹⁰; 7 × 10⁻¹⁰ | m² s⁻¹ | About half their values in water; one value for the charged components, about that of potassium lactate, so that diffusion separates no charge. |
| Stimulated flow at most; the sugar that gives half of it | 2.0; 50 | mL min⁻¹; mM | Orders of magnitude for tasting sugar. |
| Sodium and potassium in saliva, at rest and stimulated | 5 and 25; 22 and 20 | mM | Typical values; chloride is set to make each saliva neutral at its measured pH. |
| Ammonium in saliva | 2 | mM | An order of magnitude. |
| Mixing of the film during an intake | 1 | s⁻¹ | Fast mixing while sugar is in the mouth, as Dibdin (1990) assumed. |
| The sipped drink; the pocket | 100 mL of 10% sucrose over 20 min; 0.02 mol m⁻² of sugar in food particles, dissolving at 3 h⁻¹ | — | Illustrative habits of a high-sugar eater. |

---

## 9. Conditions under which reference data were obtained

Published rate constants are only meaningful alongside the conditions that
produced them. MARSE records these as metadata on each parameter set so that
a simulation states what its numbers describe.

| Field | Why it is recorded |
|---|---|
| Organism and strain | Cardinal values and $K_S$ are strain-dependent |
| Medium type (rich vs defined) | Sets which substrate is limiting and whether $K_S$ is meaningful |
| Carbon source and concentration | Determines yield and the substrate the kinetics refer to |
| Temperature | Scales every rate (theory.md §2.2) |
| pH and buffering | Scales rates; buffering determines whether pH stays constant |
| Aeration / gas phase | Sets $C_0$ for oxygen; a 5 % CO₂ atmosphere changes it |
| Culture format | Batch, chemostat, colony or flow cell — determines which model applies |
| Measurement method | Optical density, dry weight and counts are not interchangeable |

Two recurring pitfalls this metadata is designed to catch:

- **Rich versus defined media.** In rich medium the limiting resource is often
  not the one being varied — reported growth arrest can reflect exhaustion of
  catabolizable amino acids rather than of an added sugar. A Monod term
  written for the added sugar would then be modelling the wrong substrate
  entirely.
- **Optical density is not biomass.** OD is a light-scattering measurement
  whose relationship to dry mass depends on cell size and shape, both of which
  vary with growth rate. Converting OD to biomass requires a
  condition-specific calibration, recorded as its own parameter.

---

## 10. Scope

The parameters here support simulation of microbial growth and biofilm
structure under defined physical and chemical conditions, for the research
questions listed in [`docs/specification.md`](specification.md).

MARSE is a simulation framework. Its outputs are consequences of the
assumptions in [`docs/theory.md`](theory.md) and the numbers in this file.
They are not experimental measurements, and they do not support clinical,
diagnostic or treatment claims — see the scope boundaries and the explicit
list of things v1.0 does not claim in the specification.

---

## 11. Contributing a parameter

To add a value (see [`CONTRIBUTING.md`](../CONTRIBUTING.md)):

1. Give the value **with units**.
2. Cite a primary source with a DOI.
3. Record the conditions of §9.
4. State the confidence level using the **A/B/C** scale above.
5. If the value is fitted rather than measured, say so and give the data and
   procedure.

A parameter without units or a source is not accepted, regardless of how
standard it appears. "Everyone uses 0.5" is not a citation.

---

## References

Full bibliographic details for the sources cited here are listed in
[`docs/theory.md` §12](theory.md#12-references). Additional sources referenced
in this file only:

- Al-Ahmad, A., Karygianni, L., Schulze Wartenhorst, M., Bächle, M., Hellwig, E., Follo, M., Vach, K. & Han, J.-S. (2016) Bacterial adhesion and biofilm formation on yttria-stabilized, tetragonal zirconia and titanium oral implant materials with low surface roughness – an in situ study. *Journal of Medical Microbiology* **65**:596–604. [doi:10.1099/jmm.0.000267](https://doi.org/10.1099/jmm.0.000267)
- Al-Radha, A.S.D., Dymock, D., Younes, C. & O'Sullivan, D. (2012) Surface properties of titanium and zirconia dental implant materials and their effect on bacterial adhesion. *Journal of Dentistry* **40**:146–153. [doi:10.1016/j.jdent.2011.12.006](https://doi.org/10.1016/j.jdent.2011.12.006)
- Collins, L.M.C. & Dawes, C. (1987) The surface area of the adult human mouth and thickness of the salivary film covering the teeth and oral mucosa. *Journal of Dental Research* **66**:1300–1302. [doi:10.1177/00220345870660080201](https://doi.org/10.1177/00220345870660080201)
- Dawes, C., Watanabe, S., Biglow-Lecomte, P. & Dibdin, G.H. (1989) Estimation of the velocity of the salivary film at some different locations in the mouth. *Journal of Dental Research* **68**:1479–1482. [doi:10.1177/00220345890680110201](https://doi.org/10.1177/00220345890680110201)
- Oda, Y., Miura, T., Mori, G., Sasaki, H., Ito, T., Yoshinari, M. *et al.* (2020) Adhesion of streptococci to titanium and zirconia. *PLOS ONE* **15**:e0234524. [doi:10.1371/journal.pone.0234524](https://doi.org/10.1371/journal.pone.0234524)
- Presser, K.A., Ratkowsky, D.A. & Ross, T. (1997) Modelling the growth rate of *Escherichia coli* as a function of pH and lactic acid concentration. *Applied and Environmental Microbiology* **63**:2355–2360. [doi:10.1128/aem.63.6.2355-2360.1997](https://doi.org/10.1128/aem.63.6.2355-2360.1997)
- Prohaska, T., Irrgeher, J., Benefield, J. *et al.* (2022) Standard atomic weights of the elements 2021 (IUPAC Technical Report). *Pure and Applied Chemistry* **94**:573–600. [doi:10.1515/pac-2019-0603](https://doi.org/10.1515/pac-2019-0603)
- Scarano, A., Piattelli, M., Caputi, S., Favero, G.A. & Piattelli, A. (2004) Bacterial adhesion on commercially pure titanium and zirconium oxide disks: an in vivo human study. *Journal of Periodontology* **75**:292–296. [doi:10.1902/jop.2004.75.2.292](https://doi.org/10.1902/jop.2004.75.2.292)
- Sezonov, G., Joseleau-Petit, D. & D'Ari, R. (2007) *Escherichia coli* physiology in Luria-Bertani broth. *Journal of Bacteriology* **189**:8746–8749. [doi:10.1128/JB.01368-07](https://doi.org/10.1128/JB.01368-07)
- Sjollema, J., Busscher, H.J. & Weerkamp, A.H. (1988) Deposition of oral streptococci and polystyrene latices onto glass in a parallel plate flow cell. *Biofouling* **1**:101–112. [doi:10.1080/08927018809378100](https://doi.org/10.1080/08927018809378100)
- Soupene, E., van Heeswijk, W.C., Plumbridge, J. *et al.* (2003) Physiological studies of *Escherichia coli* strain MG1655: growth defects and apparent cross-regulation of gene expression. *Journal of Bacteriology* **185**:5611–5626. [doi:10.1128/JB.185.18.5611-5626.2003](https://doi.org/10.1128/JB.185.18.5611-5626.2003)
- Van Derlinden, E., Bernaerts, K. & Van Impe, J.F. (2008) Accurate estimation of cardinal growth temperatures of *Escherichia coli* from optimal dynamic experiments. *International Journal of Food Microbiology* **128**:89–100. [doi:10.1016/j.ijfoodmicro.2008.07.014](https://doi.org/10.1016/j.ijfoodmicro.2008.07.014)
