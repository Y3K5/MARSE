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

Implemented in `marse.spatial.solutes`. Both oxygen correlations refuse to
extrapolate outside their stated validity ranges (0–40 °C and 0–95 °C
respectively) rather than returning a plausible-looking wrong number.

**Corrections MARSE does not apply automatically.** Dissolved salts lower
oxygen solubility, so a culture medium holds less than pure water; a 5 % CO₂
atmosphere displaces oxygen to roughly 0.95 × the tabulated value. Both must
be stated in configuration.

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

## 7. Conditions under which reference data were obtained

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

## 8. Scope

The parameters here support simulation of microbial growth and biofilm
structure under defined physical and chemical conditions, for the research
questions listed in [`docs/specification.md`](specification.md).

MARSE is a simulation framework. Its outputs are consequences of the
assumptions in [`docs/theory.md`](theory.md) and the numbers in this file.
They are not experimental measurements, and they do not support clinical,
diagnostic or treatment claims — see the scope boundaries and the explicit
list of things v1.0 does not claim in the specification.

---

## 9. Contributing a parameter

To add a value (see [`CONTRIBUTING.md`](../CONTRIBUTING.md)):

1. Give the value **with units**.
2. Cite a primary source with a DOI.
3. Record the conditions of §7.
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

- Presser, K.A., Ratkowsky, D.A. & Ross, T. (1997) Modelling the growth rate of *Escherichia coli* as a function of pH and lactic acid concentration. *Applied and Environmental Microbiology* **63**:2355–2360. [doi:10.1128/aem.63.6.2355-2360.1997](https://doi.org/10.1128/aem.63.6.2355-2360.1997)
- Sezonov, G., Joseleau-Petit, D. & D'Ari, R. (2007) *Escherichia coli* physiology in Luria-Bertani broth. *Journal of Bacteriology* **189**:8746–8749. [doi:10.1128/JB.01368-07](https://doi.org/10.1128/JB.01368-07)
- Soupene, E., van Heeswijk, W.C., Plumbridge, J. *et al.* (2003) Physiological studies of *Escherichia coli* strain MG1655: growth defects and apparent cross-regulation of gene expression. *Journal of Bacteriology* **185**:5611–5626. [doi:10.1128/JB.185.18.5611-5626.2003](https://doi.org/10.1128/JB.185.18.5611-5626.2003)
- Van Derlinden, E., Bernaerts, K. & Van Impe, J.F. (2008) Accurate estimation of cardinal growth temperatures of *Escherichia coli* from optimal dynamic experiments. *International Journal of Food Microbiology* **128**:89–100. [doi:10.1016/j.ijfoodmicro.2008.07.014](https://doi.org/10.1016/j.ijfoodmicro.2008.07.014)
