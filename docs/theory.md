# Theory: the mathematics of MARSE

This document defines every equation MARSE solves, the assumptions behind it,
and the numerical methods used. It is the reference that
[`docs/validation.md`](validation.md) tests against and that
[`docs/parameters.md`](parameters.md) supplies numbers for.

Nothing here is specific to one organism. Model choices are provider-level
decisions recorded in each run manifest (see
[`docs/architecture.md`](architecture.md)), so a run always states which of
these equations produced it.

## 0. Conventions

### 0.1 State variables

| Symbol | Meaning | SI-consistent unit | MARSE internal unit |
|---|---|---|---|
| $X$ | biomass concentration | kg m<sup>-3</sup> | g L<sup>-1</sup> (= kg m<sup>-3</sup>) |
| $S$ | dissolved substrate concentration | mol m<sup>-3</sup> | mM (= mol m<sup>-3</sup>) |
| $C$ | dissolved oxygen concentration | mol m<sup>-3</sup> | mM |
| $\mu$ | specific growth rate | s<sup>-1</sup> | h<sup>-1</sup> |
| $t$ | time | s | h |
| $x, y, z$ | position | m | µm |
| $D$ | diffusion coefficient | m<sup>2</sup> s<sup>-1</sup> | µm<sup>2</sup> s<sup>-1</sup> |
| $T$ | temperature | K | °C |

g L<sup>-1</sup> and mM are numerically identical to kg m<sup>-3</sup> and
mol m<sup>-3</sup>, so the internal system is SI up to the length and time
scales. Length in µm and time in hours are the natural scales for a biofilm;
conversion happens once, at the configuration boundary, and the manifest
records which unit each input used.

Two mixed-unit combinations recur and are worth stating explicitly:

$$
1\ \text{µm}^2\,\text{s}^{-1} = 10^{-12}\ \text{m}^2\,\text{s}^{-1},
\qquad
1\ \text{mM} = 1\ \text{mol}\,\text{m}^{-3} = 10^{-15}\ \text{mol}\,\text{µm}^{-3}.
$$

### 0.2 Biomass units

MARSE tracks biomass as dry mass per volume, not cell count. Cell number is a
derived quantity, because dry mass per cell varies several-fold with growth
rate. Where a model needs counts, the conversion factor is an explicit,
sourced parameter, never a hidden constant.

### 0.3 What "growth rate" means

For a population of density $X$ growing without spatial structure,

$$
\frac{dX}{dt} = \mu X
\qquad\Longleftrightarrow\qquad
\mu = \frac{1}{X}\frac{dX}{dt} = \frac{d\ln X}{dt}.
$$

The doubling time follows from integrating at constant $\mu$:

$$
t_d = \frac{\ln 2}{\mu}.
$$

Implemented as `marse.microbes.growth.doubling_time` and
`specific_growth_rate`. Converting between the two is the single most common
source of factor-of-$\ln 2$ errors in the literature, so MARSE stores $\mu$
internally and converts only for display.

---

## 1. Primary models: growth under fixed conditions

Primary models describe how a population changes over time when the
environment is held constant.

### 1.1 Unrestricted (exponential) growth

$$
\frac{dX}{dt} = \mu_{\max} X
\qquad\Longrightarrow\qquad
X(t) = X_0 e^{\mu_{\max} t}.
$$

This is the limiting behaviour every other model must reproduce when
substrate is abundant and the population is far from any ceiling. It is never
valid for long: it has no mechanism to stop.

Reference solution: `marse.validation.analytical.exponential_growth`.

### 1.2 Logistic growth

The simplest closure that stops:

$$
\frac{dX}{dt} = \mu_{\max} X\left(1 - \frac{X}{X_{\max}}\right)
\qquad\Longrightarrow\qquad
X(t) = \frac{X_{\max}}{1 + \left(\dfrac{X_{\max}}{X_0} - 1\right)e^{-\mu_{\max} t}}.
$$

MARSE treats the logistic equation as a **phenomenological** curve, not a
mechanism. $X_{\max}$ is fitted, not predicted; it absorbs substrate
exhaustion, waste accumulation and space limitation into one number without
distinguishing them. It is useful as a null model and as a validation target,
and it is the wrong tool when the question is *why* growth stopped.

Reference solution: `marse.validation.analytical.logistic_growth`.

### 1.3 Monod kinetics

The mechanistic alternative: growth rate is set by the concentration of a
limiting substrate.

$$
\mu(S) = \mu_{\max}\,\frac{S}{K_S + S}.
$$

$K_S$, the half-saturation constant, is the concentration at which growth runs
at half its maximum. The form is hyperbolic: linear in $S$ when
$S \ll K_S$, saturating at $\mu_{\max}$ when $S \gg K_S$.

Two properties matter for simulation:

- **It never reaches $\mu_{\max}$.** At $S = 10 K_S$ the rate is $0.909\mu_{\max}$.
- **It is not zero at zero.** $\mu(0) = 0$ exactly, but the approach is linear,
  so tiny residual substrate still supports tiny growth. Real populations have
  a threshold substrate concentration below which they do not grow at all;
  Monod does not represent it. See §1.3.1.

Monod is an empirical fit to whole-cell behaviour that happens to share the
algebraic form of Michaelis–Menten enzyme kinetics. The resemblance is not a
derivation — $K_S$ is a population-level property that depends on transporter
expression, and it varies with growth history.

Implemented as `marse.microbes.growth.monod`. Negative concentrations, which
appear as numerical undershoot in explicit solvers, are clamped to zero rather
than producing negative growth.

### 1.3.1 The growth threshold $S_{\min}$

Subtracting maintenance and decay at rate $b$ gives a net rate that is
*negative* at low substrate and crosses zero at a finite concentration:

$$
\mu_{\text{net}}(S) = \mu_{\max}\frac{S}{K_S + S} - b,
\qquad
S_{\min} = \frac{K_S\,b}{\mu_{\max} - b},
$$

with $S_{\min} = \infty$ when $b \ge \mu_{\max}$ — the population cannot
sustain itself at any concentration. Below $S_{\min}$, substrate flux only
meets maintenance demand and there is no net growth.

This matters most in exactly the conditions MARSE is eventually aimed at.
Substrate concentrations in natural environments are orders of magnitude below
those in batch culture, and there the difference between "grows very slowly"
and "does not grow" decides which organism persists
([`modeling-landscape.md` §3.4](modeling-landscape.md#34-there-is-a-floor-the-minimum-substrate-concentration)).
Bare Monod kinetics, integrated over a long simulation, accumulates biomass
that should not exist.

Note that this is the same algebra as the chemostat break-even concentration of
§7.1, with maintenance in place of the dilution rate. Implemented as
`net_growth_rate` and `minimum_substrate_concentration`.

### 1.4 Several limiting substrates

When growth needs more than one substrate (a carbon source *and* oxygen, say),
MARSE supports two closures, chosen explicitly per run:

**Multiplicative (interactive):**

$$
\mu = \mu_{\max}\,\frac{S}{K_S + S}\cdot\frac{C}{K_C + C}.
$$

**Minimum (Liebig, non-interactive):**

$$
\mu = \mu_{\max}\,\min\!\left(\frac{S}{K_S + S},\ \frac{C}{K_C + C}\right).
$$

They differ substantially when two substrates are *both* near-limiting: the
multiplicative form compounds the two penalties, the minimum form applies only
the worse one. Neither is universally right. MARSE requires the choice to be
declared in configuration and recorded in the manifest, because a result that
flips between them is a result about the closure, not about the biology.

### 1.5 Substrate inhibition (Haldane–Andrews)

Some substrates are toxic in excess:

$$
\mu(S) = \mu_{\max}\,\frac{S}{K_S + S + \dfrac{S^2}{K_I}}.
$$

The rate is non-monotonic, peaking at

$$
S^{*} = \sqrt{K_S K_I},
\qquad
\mu(S^{*}) = \frac{\mu_{\max}}{1 + 2\sqrt{K_S/K_I}}.
$$

Note that the peak is strictly below $\mu_{\max}$ — the nominal $\mu_{\max}$ in
this model is not attainable, which makes parameters fitted with Haldane
non-comparable with parameters fitted with Monod. As $K_I \to \infty$ the model
reduces to Monod exactly.

This non-monotonicity has a real consequence: in a chemostat it creates two
steady states at the same dilution rate, one stable and one unstable, so the
outcome depends on the initial condition (Andrews 1968).

Implemented as `marse.microbes.growth.haldane`.

### 1.6 Lag phase: the Baranyi–Roberts model

Cells transferred into new conditions do not grow immediately. Baranyi &
Roberts (1994) model this with an *adjustment function* representing the
physiological state of the inoculum. Writing $y = \ln X$:

$$
\frac{dy}{dt} = \mu_{\max}\,\underbrace{\frac{Q(t)}{1 + Q(t)}}_{\text{lag}}\,
\underbrace{\left[1 - e^{m(y - y_{\max})}\right]}_{\text{braking}},
\qquad
\frac{dQ}{dt} = \mu_{\max} Q .
$$

$Q$ is the concentration of a hypothetical limiting "work to be done" that
must accumulate before growth proceeds. Under constant conditions this
integrates to a closed form:

$$
y(t) = y_0 + \mu_{\max} A(t)
- \frac{1}{m}\ln\!\left[1 + \frac{e^{m\mu_{\max}A(t)} - 1}{e^{m(y_{\max} - y_0)}}\right],
$$

$$
A(t) = t + \frac{1}{\mu_{\max}}
\ln\!\left(e^{-\mu_{\max} t} + e^{-h_0} - e^{-\mu_{\max} t - h_0}\right),
\qquad
h_0 = \mu_{\max}\lambda .
$$

Here $\lambda$ is the lag duration and $m$ shapes the approach to stationary
phase ($m = 1$ gives the logistic braking term).

Three properties make this the preferred primary model in MARSE:

1. **$h_0$ is the physiologically meaningful quantity**, not $\lambda$. It is
   dimensionless "work" carried by the inoculum, and it is approximately
   invariant when the same inoculum is shifted into different temperatures —
   which means lag under changing conditions is predictable rather than
   refitted.
2. The **differential form extends to time-varying conditions**, so the same
   model handles perturbation experiments. The closed form above applies only
   when conditions are constant.
3. It reduces cleanly: $A(t) \to t$ as $\lambda \to 0$, and
   $y \to y_0 + \mu_{\max}t$ far from $y_{\max}$.

Asymptotically the exponential phase is the no-lag line shifted right by
exactly $\lambda$, which is the behaviour the test suite checks.

Implemented as `marse.microbes.growth.baranyi_roberts`, evaluated with
`logaddexp`/`log1p` so that the exponentials do not overflow for large
$m\mu_{\max}A$.

### 1.7 What MARSE does not model at this layer

Death, lysis, viable-but-nonculturable states and sporulation are separate
processes with their own rates. Folding them into a "net growth rate" hides
the distinction between a population that is not growing and one that is
growing and dying at equal rates — which behave identically in the model and
completely differently under perturbation. MARSE keeps them separate (Phase 4).

---

## 2. Secondary models: how the environment scales growth

Secondary models supply $\mu_{\max}$ as a function of temperature, pH and
other conditions.

### 2.1 The gamma concept

Zwietering et al. (1992) proposed that environmental factors act
independently and multiplicatively on the optimal rate:

$$
\mu_{\max}(T, \mathrm{pH}, a_w, \ldots)
= \mu_{\mathrm{opt}}\cdot\gamma_T(T)\cdot\gamma_{\mathrm{pH}}(\mathrm{pH})\cdot\gamma_{a_w}(a_w)\cdots
$$

with each $\gamma \in [0, 1]$ and $\gamma = 1$ at the optimum. This is
convenient — each factor is calibrated from separate single-factor
experiments — and it is an approximation.

**It is known to fail near the growth limits.** Baka et al. (2013) showed for
*E. coli* K-12 that the cardinal temperatures themselves shift with pH,
following a parabolic trend, which the gamma hypothesis forbids. The
independence assumption is reasonable in the interior of the growth region and
degrades as any factor approaches its limit.

MARSE therefore:

- uses the gamma form as the default, because it is calibratable;
- records in the manifest that it was used;
- states this limitation in the output metadata rather than in a footnote;
- leaves the provider interface open for interaction models.

Conclusions drawn near a growth boundary should not rest on the gamma form.

### 2.2 Cardinal temperature model with inflection (CTMI)

Rosso et al. (1993). For $T_{\min} < T < T_{\max}$:

$$
\gamma_T(T) =
\frac{(T - T_{\max})(T - T_{\min})^2}
{(T_{\mathrm{opt}} - T_{\min})\Big[(T_{\mathrm{opt}} - T_{\min})(T - T_{\mathrm{opt}})
- (T_{\mathrm{opt}} - T_{\max})(T_{\mathrm{opt}} + T_{\min} - 2T)\Big]},
$$

and $\gamma_T = 0$ outside. The three parameters are directly interpretable —
they are the temperatures at which growth stops and is fastest — which is why
MARSE prefers this over models whose parameters are fitting artefacts.

Two structural facts:

- Because only temperature *differences* appear, °C and K give identical
  results. The test suite asserts this.
- The formula is well posed only when
  $T_{\mathrm{opt}} \ge \tfrac{1}{2}(T_{\min} + T_{\max})$, i.e. the optimum
  lies nearer the maximum than the minimum. This holds for essentially all
  measured microbial data — the fall-off above the optimum is steep, the rise
  below it is gradual — and MARSE rejects parameter sets that violate it
  rather than silently producing a curve with $\gamma > 1$.

The asymmetry is not a modelling convenience; it reflects that the high-
temperature limit is set by protein denaturation, which is cooperative and
therefore abrupt, while the low-temperature limit is set by reaction rates
slowing smoothly.

Implemented as `marse.microbes.cardinal.cardinal_temperature`.

### 2.3 Cardinal pH model (CPM)

Rosso et al. (1995), same philosophy:

$$
\gamma_{\mathrm{pH}}(\mathrm{pH}) =
\frac{(\mathrm{pH} - \mathrm{pH}_{\min})(\mathrm{pH} - \mathrm{pH}_{\max})}
{(\mathrm{pH} - \mathrm{pH}_{\min})(\mathrm{pH} - \mathrm{pH}_{\max}) - (\mathrm{pH} - \mathrm{pH}_{\mathrm{opt}})^2},
$$

zero outside $[\mathrm{pH}_{\min}, \mathrm{pH}_{\max}]$. Symmetric about
$\mathrm{pH}_{\mathrm{opt}}$, unlike the temperature model.

A caution that matters for biofilms specifically: pH is a *local* variable. In
a metabolically active biofilm producing organic acids, local pH can differ
from bulk pH by more than a unit. Applying a bulk-pH $\gamma$ factor uniformly
across a biofilm is an approximation that MARSE makes explicit; modelling the
acid as a diffusing solute with local pH is the more defensible route once
Phase 3 lands.

Implemented as `marse.microbes.cardinal.cardinal_ph`.

### 2.4 Ratkowsky square-root model

Ratkowsky et al. (1982, 1983) give an alternative:

$$
\sqrt{\mu} = b\,(T - T_{\min})\left[1 - e^{c(T - T_{\max})}\right].
$$

Below the optimum the exponential term is negligible and this reduces to the
linear relation $\sqrt{\mu} = b(T - T_{\min})$, which is the form most often
plotted in the food-microbiology literature. Note it returns the rate itself,
not a $\gamma$ factor — $b$ carries units.

The Arrhenius equation, by contrast, does *not* describe microbial growth over
the biokinetic range: plots of $\ln\mu$ against $1/T$ curve rather than
straightening (Ratkowsky et al. 1982). MARSE does not offer an Arrhenius
growth provider.

Implemented as `marse.microbes.cardinal.ratkowsky`.

---

## 3. Stoichiometry: linking growth to consumption

### 3.1 Yield

The growth yield relates biomass produced to substrate consumed:

$$
Y_{X/S} = \frac{\Delta X}{-\Delta S}.
$$

Under pure growth with no maintenance, substrate uptake is

$$
q_S = \frac{\mu}{Y_{X/S}}
\qquad [\text{substrate} \cdot \text{biomass}^{-1} \cdot \text{time}^{-1}].
$$

### 3.2 Maintenance (Pirt)

Cells consume substrate even when not growing — to maintain gradients, turn
over protein, and stay alive. Pirt (1965):

$$
q_S = \frac{\mu}{Y_{X/S}^{\mathrm{true}}} + m_S .
$$

This has an immediate and important consequence. The *observed* yield is

$$
\frac{1}{Y^{\mathrm{obs}}} = \frac{1}{Y^{\mathrm{true}}} + \frac{m_S}{\mu},
$$

so **observed yield falls as growth slows**. A yield measured in a fast batch
culture systematically over-predicts biomass in a slow-growing biofilm, where
much of the population grows slowly and maintenance dominates. Using a
batch-measured yield in a biofilm simulation is a specific, recurring error
that this term exists to prevent.

The volumetric consumption rate is then

$$
r_S = \left(\frac{\mu}{Y^{\mathrm{true}}} + m_S\right) X .
$$

Implemented as `marse.microbes.growth.substrate_uptake_rate`.

### 3.3 Product formation (Luedeking–Piret)

Metabolites, extracellular polymeric substances and acids are produced with
both growth-associated and non-growth-associated components:

$$
r_P = \left(\alpha\mu + \beta\right)X .
$$

$\alpha$ is growth-associated (product made in proportion to new biomass),
$\beta$ is not (product made by existing biomass regardless). MARSE uses this
form for matrix production and for cross-feeding metabolites, where the two
components behave very differently: a $\beta$-dominated product keeps
accumulating in a stationary biofilm, an $\alpha$-dominated one does not.

Implemented as `marse.microbes.growth.product_formation_rate`.

### 3.4 Mass balance as an invariant

In a closed batch system with no maintenance, the combination

$$
X + Y_{X/S}\,S = X_0 + Y_{X/S}\,S_0
$$

is conserved exactly. MARSE's test suite asserts this to machine precision
during integration, because a solver that violates it is producing biomass
from nothing — the single most damaging class of silent error in this kind of
model. The final biomass of a batch culture follows immediately:

$$
X_\infty = X_0 + Y_{X/S} S_0 .
$$

Implemented as `marse.validation.analytical.batch_final_biomass`.

### 3.5 Integrated Monod batch solution

Combining Monod growth with the mass balance yields an exact implicit
solution for batch culture — the time at which substrate has fallen to $S$:

$$
\mu_{\max}t = \frac{K_S}{C}\ln\!\frac{S_0}{S}
+ \left(1 + \frac{K_S}{C}\right)\ln\!\frac{X}{X_0},
\qquad
C = S_0 + \frac{X_0}{Y},
\quad
X = X_0 + Y(S_0 - S).
$$

This is a genuine analytical reference for a nonlinear coupled system, and it
is the strongest correctness test available for the growth–uptake coupling.
The test suite verifies it against fourth-order Runge–Kutta integration to a
relative tolerance of $10^{-7}$.

Implemented as `marse.validation.analytical.monod_batch_time`.

### 3.6 Composition, continuity and the degree of reduction

Sections 3.1–3.5 conserve mass for one substrate and one biomass. A reaction
network (configuration schema version 2, [`docs/networks.md`](networks.md))
makes conservation a property of the configuration itself, for any number of
components and processes. It follows the continuity check of the IWA
activated-sludge and biofilm models (Henze et al. 2000).

Each component $j$ has a chemical formula
$\mathrm{C}_{c_j}\mathrm{H}_{h_j}\mathrm{O}_{o_j}\mathrm{N}_{n_j}$ with charge
$z_j$, counted per mol of that formula unit (per C-mol for biomass written per
carbon atom). Its composition in the three conserved quantities is

$$
I_{j,\mathrm{C}} = c_j,
\qquad
I_{j,\mathrm{N}} = n_j,
\qquad
I_{j,e} = \gamma_j = 4c_j + h_j - 2o_j - 3n_j - z_j .
$$

$\gamma_j$ is the **degree of reduction** (Roels 1983): the electrons released
when the component is oxidised completely to CO₂, H₂O, NH₄⁺ and H⁺, the
reference state in which it is zero. Glucose has 24, biomass
CH<sub>1.8</sub>O<sub>0.5</sub>N<sub>0.2</sub> has 4.2 (Heijnen and van Dijken
1992), and O₂, which accepts electrons, has −4. Protonation leaves it
unchanged: acetic acid and acetate both have 8. Because one mol of electrons
reduces a quarter mol of O₂, $8\gamma_j$ is the chemical oxygen demand in
grams per mol. COD, the unit of the activated-sludge models, is therefore the
electron balance written in mass.

A process $p$ is a row $\nu_{pj}$ of the stoichiometric (Gujer) matrix,
negative for what it consumes and positive for what it produces. It conserves
matter if and only if

$$
\sum_j \nu_{pj}\, I_{jk} = 0 \qquad \text{for } k = \mathrm{C},\ \mathrm{N},\ e .
$$

MARSE refuses at load time any process that fails this, naming the process
and the quantity.

**Water and protons close the rest.** Oxygen, hydrogen and charge are not
tracked. Define the implicit water and proton coefficients from the oxygen
and charge balances,

$$
\nu_{\mathrm{H_2O}} = -\sum_j \nu_{pj}\, o_j,
\qquad
\nu_{\mathrm{H^+}} = -\sum_j \nu_{pj}\, z_j .
$$

The hydrogen balance then holds automatically. The three conditions give
$\sum_j \nu_{pj}(4c_j + h_j - 2o_j - 3n_j - z_j) = 0$ with
$\sum_j \nu_{pj} c_j = \sum_j \nu_{pj} n_j = 0$, so

$$
\sum_j \nu_{pj}\, h_j = \sum_j \nu_{pj}\,(2o_j + z_j) = -2\nu_{\mathrm{H_2O}} - \nu_{\mathrm{H^+}} .
$$

That is exactly the hydrogen balance with water and protons included. So every
row that satisfies the three conditions is a complete chemical equation once
H₂O and H⁺ are added, and no row that fails them can be completed. This is why
three quantities suffice, and why water and protons must never be listed as
components: they hold no carbon, nitrogen or electrons, so no balance would
constrain them.

**Balancing.** A process gives some coefficients and names a set $B$ of at most
three components whose coefficients $x_b$ the balances determine:

$$
\sum_{b \in B} x_b\, I_{bk} = -\sum_{j \notin B} \nu_{pj}\, I_{jk}
\qquad \text{for } k = \mathrm{C},\ \mathrm{N},\ e .
$$

MARSE requires exactly one solution and finds it by Gauss–Jordan elimination
in rational arithmetic. For growth of biomass $X$ on substrate $S$ with yield
$Y_{X/S}$ (mol per mol), the row is written per mol of biomass formed:
$\nu_X = 1$, $\nu_S = -1/Y_{X/S}$, and $\nu_P = Y_{P/S}/Y_{X/S}$ for each
product with a stated yield. The electron acceptor, carbon dioxide and the
nitrogen source are usually in $B$. This is the electron-balance calculation
of the textbooks; the test suite checks it against the half-reaction method of
Rittmann and McCarty (2001, ch. 2).

**Exact, then rounded once.** Coefficients are read as the decimals written in
the configuration, so the check is exact, not approximate: a coefficient
copied as 5.999 instead of 6 is refused. The engine receives each coefficient
rounded once to double precision, so $\sum_j \nu_{pj} I_{jk}$ is zero to within
a few units in the last place. The runtime ledger of §9.5 is built on that
margin.

Implemented as `marse.schemas.network`.

---

## 4. Transport

### 4.1 Diffusion

Solutes move by Fick's second law, with a reaction term:

$$
\frac{\partial C}{\partial t} = \nabla\cdot\left(D\nabla C\right) - r(C, X).
$$

For constant $D$ this is $\partial_t C = D\nabla^2 C - r$.

The diffusive time scale over a distance $L$ is

$$
\tau_D \sim \frac{L^2}{D}.
$$

The quadratic dependence is the single most important scaling fact in biofilm
modelling: doubling the thickness quadruples the transport time.

### 4.2 Point-source reference solution

For an instantaneous release of amount $M$ at the origin in two dimensions
with no reaction:

$$
C(r, t) = \frac{M}{4\pi D t}\,e^{-r^2/(4Dt)} .
$$

This is the analytical reference for verifying a diffusion solver
independently of any biology — it conserves mass exactly and satisfies the
diffusion equation exactly, so a solver that reproduces it is transporting
correctly.

Implemented as `marse.validation.analytical.point_source_diffusion_2d`.

### 4.3 Effective diffusivity in biofilms

Diffusion inside a biofilm is slower than in water. MARSE parameterises this
by a relative effective diffusivity:

$$
D_e = f\cdot D_{\mathrm{aq}}, \qquad 0 < f \le 1 .
$$

Stewart (1998) reviewed measured values and found $f$ depends mainly on solute
physical chemistry rather than on the organism, reporting mean relative
effective diffusive permeabilities of approximately 0.56 for inorganic
ions, 0.43 for small nonpolar solutes (molecular weight ≤ 44, which includes
oxygen at 32), and 0.29 for larger organic solutes.

$f$ is a **lumped parameter**, not a material constant: it absorbs matrix
density, porosity, cell packing and binding. It varies within a single biofilm
and over its lifetime. MARSE treats it as a per-solute, per-run parameter with
a recorded source, and sensitivity to $f$ should be reported for any
conclusion that depends on gradients.

### 4.4 Oxygen: solubility

Oxygen concentration at an air interface is set by solubility, which falls
markedly with temperature. MARSE uses the Benson & Krause (1984) correlation
for fresh water in equilibrium with air at 1 atm:

$$
\ln C^{*} = -139.34411 + \frac{1.575701\times10^{5}}{T_K}
- \frac{6.642308\times10^{7}}{T_K^{2}}
+ \frac{1.243800\times10^{10}}{T_K^{3}}
- \frac{8.621949\times10^{11}}{T_K^{4}},
$$

with $C^{*}$ in mg L<sup>-1</sup> and $T_K$ in kelvin; valid 0–40 °C.

| $T$ (°C) | $C^{*}$ (mg L<sup>-1</sup>) | $C^{*}$ (mM) |
|---|---|---|
| 20 | 9.09 | 0.284 |
| 25 | 8.26 | 0.258 |
| 30 | 7.56 | 0.236 |
| 37 | 6.73 | 0.210 |

Two corrections MARSE requires the user to make explicit rather than applying
silently: **dissolved salts lower solubility** (culture media hold measurably
less oxygen than pure water), and a **5 % CO₂ atmosphere** displaces oxygen
proportionally, giving roughly $0.95\,C^{*}$.

Implemented as `marse.spatial.solutes.oxygen_saturation_mg_per_l`.

### 4.5 Oxygen: diffusivity

Han & Bartels (1996), measured 0–95 °C:

$$
\log_{10}\!\left(\frac{D}{\mathrm{cm^2\,s^{-1}}}\right)
= -4.410 + \frac{773.8}{T_K} - \left(\frac{506.4}{T_K}\right)^{2}.
$$

At 25 °C this gives $2.00\times10^{-9}$ m² s<sup>-1</sup> (2000 µm² s<sup>-1</sup>);
at 37 °C, $2.62\times10^{-9}$ m² s<sup>-1</sup> (2624 µm² s<sup>-1</sup>).

Note the deliberate asymmetry with §4.4: **warming water increases the
diffusivity of oxygen while decreasing how much of it dissolves.** The two
effects act in opposite directions on delivery, so neither can be neglected
when comparing runs at different temperatures.

Implemented as `marse.spatial.solutes.oxygen_diffusivity_m2_per_s`.

### 4.6 Temperature scaling for other solutes

Where a diffusivity is known at one temperature only, MARSE scales it with the
Stokes–Einstein relation,

$$
D \propto \frac{T}{\eta(T)},
\qquad\text{so}\qquad
D(T_2) = D(T_1)\cdot\frac{T_2}{T_1}\cdot\frac{\eta(T_1)}{\eta(T_2)},
$$

using the IAPWS 2008 formulation for water viscosity $\eta$ (Huber et al.
2009). For glucose, taking 600 µm² s<sup>-1</sup> at 25 °C gives
approximately 800 µm² s<sup>-1</sup> at 37 °C.

This scaling assumes a rigid sphere in a continuum solvent. It is adequate for
small solutes over a 10–20 K span and is flagged as an approximation in the
manifest whenever it is applied.

---

## 5. Reaction–diffusion coupling

### 5.1 Penetration depth

The central quantitative question in biofilm modelling: how deep does a solute
get before it is consumed?

For **zero-order** uptake (rate $k_0$, independent of concentration — the
relevant limit when $C \gg K_C$ near the surface) in a slab with surface
concentration $C_0$, the steady state $D_e C'' = k_0$ with $C(0) = C_0$ and
$C(\delta) = C'(\delta) = 0$ gives a parabolic profile reaching zero at

$$
\boxed{\ \delta = \sqrt{\frac{2 D_e C_0}{k_0}}\ }
\qquad\text{with}\qquad
C(z) = C_0\left(1 - \frac{z}{\delta}\right)^2 .
$$

Beyond $\delta$ the solute is absent and that region is anoxic (or starved).

This single expression explains stratification: because $\delta$ scales as the
*square root* of supply and inversely as the square root of demand, a biofilm
that is metabolically active is necessarily stratified once it is thicker than
a few tens of micrometres. It is consistent with the oxygen penetration
depths of tens of micrometres reported for dense colony biofilms (Walters
et al. 2003; Werner et al. 2004), against total thicknesses several times
larger.

Implemented as `marse.validation.analytical.zero_order_penetration_depth`,
and verified in the test suite against a Newton-solved nonlinear
reaction–diffusion system with Monod uptake at small $K_C$ — including a check
that the surface flux equals $k_0\delta$, i.e. that all uptake is fed by
transport across the surface.

### 5.2 Worked example

Oxygen in a biofilm at 37 °C, with $f = 0.43$ (§4.3):

$$
D_e = 0.43 \times 2624 = 1128\ \text{µm}^2\,\text{s}^{-1},
\qquad
C_0 = 0.210\ \text{mM}.
$$

For a penetration depth of 50 µm, the required volumetric uptake is

$$
k_0 = \frac{2 D_e C_0}{\delta^2} = 0.19\ \text{mM s}^{-1} \approx 683\ \text{mM h}^{-1},
$$

which at a specific oxygen uptake rate of order 20–27 mmol gDW<sup>-1</sup> h<sup>-1</sup>
corresponds to roughly 25 g L<sup>-1</sup> dry biomass — a plausible density
for a dense biofilm. The parameters are mutually consistent, which is the kind
of order-of-magnitude check MARSE expects before a configuration is run.

The corresponding diffusive equilibration time across 100 µm is
$\tau_D = 100^2/1128 \approx 8.9$ s.

### 5.3 Time-scale separation

Compare that 8.9 s with a doubling time of 20–60 minutes. Solute transport
equilibrates **two to three orders of magnitude faster** than biomass changes.

This justifies the standard **quasi-steady-state approximation**: hold biomass
fixed, solve $\nabla\cdot(D_e\nabla C) = r(C, X)$ to steady state, then advance
biomass over a much longer time step using the resulting field. This is what
makes multi-day biofilm simulations tractable, and it is the approach taken by
established biofilm models (Wanner & Gujer 1986; Picioreanu et al. 1998;
Kreft et al. 2001; Lardon et al. 2011).

MARSE will make the separation **checkable rather than assumed**: the ratio
$\tau_D/\tau_{\text{growth}}$ is computed and recorded, and a run whose ratio
exceeds a configured threshold is flagged in its manifest. The approximation
fails during fast transients — precisely the perturbation experiments MARSE is
built for — so it cannot be silently trusted.

---

## 6. Spatial biomass dynamics

### 6.1 Biomass balance

For biomass density $X_i$ of species $i$:

$$
\frac{\partial X_i}{\partial t}
= \underbrace{\mu_i X_i}_{\text{growth}}
- \underbrace{b_i X_i}_{\text{decay}}
- \underbrace{\nabla\cdot(\mathbf{u}X_i)}_{\text{spreading}}
+ \underbrace{\Sigma_i}_{\text{attachment, detachment}} .
$$

Biomass does not diffuse. It is displaced: as cells divide in a confined
space, the accumulating volume pushes neighbouring biomass outward. The
advective velocity $\mathbf{u}$ is determined by a constraint — that local
biomass density cannot exceed a maximum packing — rather than by a
constitutive law. The three standard resolutions are continuum displacement
(Wanner & Gujer), discrete cellular-automaton shoving (Picioreanu et al.), and
individual-based mechanical relaxation (Kreft et al.). Phase 0 has not fixed
this choice; it is the main open question in
[`docs/specification.md`](specification.md).

### 6.2 Colony radial expansion

A surface colony fed from its edge grows radially at a constant rate once
established (Pirt 1967):

$$
\frac{dR}{dt} = k_r \quad\Longrightarrow\quad R(t) = R_0 + k_r t .
$$

The mechanism is direct: only an annulus of width $w$ at the periphery has
access to nutrient, so the growing mass is proportional to the perimeter,
$dA/dt \propto 2\pi R w$, and with $A = \pi R^2$ this gives constant $dR/dt$.

This makes radial expansion an unusually good validation target. **Linear
radius growth is a signature of diffusion-limited peripheral growth**, and a
simulation that produces exponential radial expansion has the transport
coupling wrong — the error is visible in the shape of the curve, not just its
magnitude, and no parameter adjustment can disguise it.

### 6.3 Detachment

Biomass leaves the biofilm by erosion (continuous, surface-proportional) and
sloughing (stochastic, large). A common erosion closure makes the rate
depend on thickness $L_f$:

$$
r_{\mathrm{det}} = k_{\mathrm{det}} L_f^{2},
$$

which produces a steady-state thickness when balanced against growth. MARSE
treats the exponent as a declared model choice, because the resulting
steady-state thickness is sensitive to it and it is poorly constrained by
data.

---

## 7. Competition and coexistence

### 7.1 The chemostat and break-even concentration

In a well-mixed chemostat at dilution rate $D$, a species persists only if its
growth balances washout, $\mu(S) = D$. For Monod kinetics this defines the
**break-even substrate concentration**:

$$
\lambda_i = \frac{K_{S,i}D}{\mu_{\max,i} - D},
\qquad
\lambda_i = \infty \text{ if } \mu_{\max,i} \le D .
$$

### 7.2 The $R^{*}$ rule

With a single limiting substrate, the species with the **lowest** $\lambda$
drives the substrate down to its own break-even level and excludes all others
(Hsu, Hubbell & Waltman 1977). The steady state is

$$
S^{*} = \lambda_{\min},
\qquad
X^{*} = Y\left(S_{\mathrm{in}} - \lambda_{\min}\right).
$$

This is a genuinely predictive, falsifiable result, and it has a
counter-intuitive consequence MARSE's test suite exercises directly: **the
winner depends on the dilution rate.** A species with low $\mu_{\max}$ but
very low $K_S$ (an efficient scavenger) wins at low dilution rates, while a
species with high $\mu_{\max}$ but high $K_S$ wins at high dilution rates. The
test runs the full three-species ODE system at two dilution rates and confirms
that the winner switches, that the residual substrate matches
$\lambda_{\min}$, and that the survivor's biomass matches $Y(S_{\mathrm{in}} - \lambda_{\min})$.

Implemented as `marse.validation.analytical.chemostat_break_even`.

### 7.3 Why spatial structure changes this

Competitive exclusion is a well-mixed result. In a biofilm, spatial
segregation, cross-feeding and gradient partitioning permit coexistence that
the chemostat forbids. This contrast — exclusion in the mixed case,
coexistence in the spatial case, *with the same kinetic parameters* — is one
of the clearest demonstrations that spatial structure is doing scientific
work, and it is the target of validation cases V5 and V6.

---

## 8. Adaptation as explicit state transitions

MARSE represents phenotypic adaptation as transitions between a finite,
declared set of states (planktonic, attached, stressed, slow-growing,
tolerant), each with its own parameters:

$$
\frac{dX_a}{dt} = \mu_a X_a + \sum_{b \ne a}\left(k_{b\to a}X_b - k_{a\to b}X_a\right).
$$

Transition rates may depend on local conditions, $k_{a\to b} = k_{a\to b}(\mathbf{c})$.

This is **not** a model of genomic evolution, and MARSE does not claim it is.
It is a model of reversible physiological switching, of the kind demonstrated
for bacterial persistence, where a subpopulation switches between normally
growing and slow-growing states at measurable rates (Balaban et al. 2004).

The design constraints are deliberate: the state set is finite and declared in
configuration; every transition rate has units, a source and a confidence
note; and switching is auditable, so a simulation can always report what
fraction of the population is in which state and why it moved. A model that
can invent new states to fit data is not falsifiable.

---

## 9. Numerical methods

### 9.1 Explicit diffusion and its stability limit

The forward-Euler discretisation of $\partial_t C = D\nabla^2 C$ on a uniform
grid of spacing $h$ in $d$ dimensions is stable only when

$$
\Delta t \le \frac{h^2}{2 d D}.
$$

This is restrictive at biofilm resolution. For oxygen in water at 37 °C
($D = 2624$ µm² s<sup>-1</sup>) on a 5 µm grid in two dimensions:

$$
\Delta t \le \frac{5^2}{4 \times 2624} \approx 2.4\ \text{ms}.
$$

Milliseconds, for a process to be simulated over days. Explicit stepping of
the full coupled system is therefore not viable, which is exactly why the
quasi-steady-state approach of §5.3 is used: the solute field is solved to
steady state directly, and only biomass is marched in time.

MARSE computes this bound from the run's own parameters and refuses a
configuration that violates it, rather than producing plausible-looking
oscillating garbage.

### 9.2 Steady-state solution

The steady solute field solves a nonlinear system (Monod uptake makes it
nonlinear in $C$). MARSE's reference implementation uses Newton iteration with
a tridiagonal (Thomas-algorithm) solve in one dimension — the approach the
test suite uses to verify §5.1 — extending to iterative solvers in higher
dimensions.

### 9.3 Operator splitting and order of updates

Within a time step, providers are applied in a **fixed order declared in
configuration and recorded in the manifest**. Splitting introduces an error of
$O(\Delta t)$ for sequential (Lie) splitting, and — critically — a different
order gives a different answer. Leaving the order implicit makes a run
irreproducible even when every parameter is identical, which is why MARSE
treats it as part of the recorded state rather than an implementation detail.

### 9.4 Randomness

All stochasticity comes from seeded generators owned by the core. Each
provider receives an independent stream derived from the run seed, so adding a
provider does not perturb the random sequence seen by the others — without
this, adding a diagnostic changes the trajectory. The seed is recorded in the
manifest and a replay reproduces the run exactly, or within stated tolerances
for ensembles.

### 9.5 Conservation checks as runtime invariants

Mass balance (§3.4) is checked during integration, not only in tests. A run
that loses or creates mass beyond tolerance is reported as failed rather than
returned as a result.

---

## 10. Assumptions and limitations

Stated plainly, because a simulation's credibility rests on what it admits it
cannot do:

1. **Continuum biomass.** Densities, not individuals. Invalid when a few cells
   determine the outcome (early attachment, rare switching events).
2. **Gamma independence.** Temperature and pH act independently — known to
   fail near growth limits (§2.1).
3. **Quasi-steady solutes.** Valid given the time-scale separation of §5.3,
   and invalid during fast transients, which is when perturbation experiments
   are most interesting.
4. **Lumped effective diffusivity.** One $f$ per solute for a structure that
   is genuinely heterogeneous (§4.3).
5. **Fixed phenotype set.** No novel states emerge (§8).
6. **Parameters from other conditions.** A $K_S$ measured in a chemostat on
   one strain in one medium is being applied elsewhere. This is the dominant
   uncertainty in practice, ahead of any numerical concern — $K_S$ has been
   reported to shift ~170-fold in one strain on one substrate through
   physiological adaptation alone
   ([`modeling-landscape.md` §3.2](modeling-landscape.md#32-the-half-saturation-constant-is-not-a-property-of-a-species)).
7. **Laboratory conditions, not natural ones.** Single limiting substrate,
   abundant nutrients, fast growth. Natural communities grow on many substrates
   at once, orders of magnitude more slowly, with large dormant fractions. v1.0
   does not represent mixed-substrate kinetics or dormancy, so it is not a model
   of a soil or marine community ([`modeling-landscape.md` §3](modeling-landscape.md#3-natural-states-versus-laboratory-states)).
8. **No host.** No immune cells, no tissue, no clinical inference (see
   [`docs/specification.md`](specification.md)).
9. **Carbon, nitrogen and electrons only.** Reaction networks (§3.6) balance
   these three; sulphur, phosphorus and metals are not yet balanced, and
   formulas containing them are refused rather than half-checked.

**Simulation output is not experimental evidence.** MARSE produces
consequences of stated assumptions. Conclusions are phrased as "under these
modelled conditions and assumptions", and the manifest exists so that any
reader can see exactly what those were.

---

## 11. Mapping to validation cases

| Case | Tests | Reference |
|---|---|---|
| V1 | Growth law, carrying capacity | §1.1, §1.2 |
| V2 | Resource limitation, yield, mass balance | §1.3, §3.4, §3.5 |
| V3 | Diffusion solver, penetration depth | §4.2, §5.1 |
| V4 | Attachment, biofilm initiation | §6.1, §6.3 |
| V5 | Two-species competition, $R^{*}$ rule | §7.1, §7.2 |
| V6 | Cross-feeding, coexistence | §3.3, §7.3 |
| V7 | Perturbation and recovery | §2, §8 |
| V8 | Reproducibility from manifest | §9.3, §9.4 |

Full definitions in [`docs/validation.md`](validation.md).

---

## 12. References

Numerical values taken from these sources, with confidence notes, are
tabulated in [`docs/parameters.md`](parameters.md).

1. Andrews, J.F. (1968) A mathematical model for the continuous culture of microorganisms utilizing inhibitory substrates. *Biotechnology and Bioengineering* **10**:707–723. [doi:10.1002/bit.260100602](https://doi.org/10.1002/bit.260100602)
2. Baka, M., Van Derlinden, E., Boons, K., Mertens, L. & Van Impe, J.F. (2013) Impact of pH on the cardinal temperatures of *E. coli* K12: evaluation of the gamma hypothesis. *Food Control* **29**:328–335. [doi:10.1016/j.foodcont.2012.04.022](https://doi.org/10.1016/j.foodcont.2012.04.022)
3. Balaban, N.Q., Merrin, J., Chait, R., Kowalik, L. & Leibler, S. (2004) Bacterial persistence as a phenotypic switch. *Science* **305**:1622–1625. [doi:10.1126/science.1099390](https://doi.org/10.1126/science.1099390)
4. Baranyi, J. & Roberts, T.A. (1994) A dynamic approach to predicting bacterial growth in food. *International Journal of Food Microbiology* **23**:277–294. [doi:10.1016/0168-1605(94)90157-0](https://doi.org/10.1016/0168-1605(94)90157-0)
5. Benson, B.B. & Krause, D. (1984) The concentration and isotopic fractionation of oxygen dissolved in freshwater and seawater in equilibrium with the atmosphere. *Limnology and Oceanography* **29**:620–632. [doi:10.4319/lo.1984.29.3.0620](https://doi.org/10.4319/lo.1984.29.3.0620)
6. Han, P. & Bartels, D.M. (1996) Temperature dependence of oxygen diffusion in H₂O and D₂O. *Journal of Physical Chemistry* **100**:5597–5602. [doi:10.1021/jp952903y](https://doi.org/10.1021/jp952903y)
7. Heijnen, J.J. & van Dijken, J.P. (1992) In search of a thermodynamic description of biomass yields for the chemotrophic growth of microorganisms. *Biotechnology and Bioengineering* **39**:833–858. [doi:10.1002/bit.260390806](https://doi.org/10.1002/bit.260390806)
8. Henze, M., Gujer, W., Mino, T. & van Loosdrecht, M.C.M. (2000) *Activated Sludge Models ASM1, ASM2, ASM2d and ASM3.* IWA Scientific and Technical Report No. 9. IWA Publishing, London.
9. Hsu, S.-B., Hubbell, S.P. & Waltman, P. (1977) A mathematical theory for single-nutrient competition in continuous cultures of micro-organisms. *SIAM Journal on Applied Mathematics* **32**:366–383. [doi:10.1137/0132030](https://doi.org/10.1137/0132030)
10. Huber, M.L., Perkins, R.A., Laesecke, A. *et al.* (2009) New international formulation for the viscosity of H₂O. *Journal of Physical and Chemical Reference Data* **38**:101–125. [doi:10.1063/1.3088050](https://doi.org/10.1063/1.3088050)
11. Kovárová-Kovar, K. & Egli, T. (1998) Growth kinetics of suspended microbial cells: from single-substrate-controlled growth to mixed-substrate kinetics. *Microbiology and Molecular Biology Reviews* **62**:646–666. [doi:10.1128/mmbr.62.3.646-666.1998](https://doi.org/10.1128/mmbr.62.3.646-666.1998)
12. Kreft, J.-U., Picioreanu, C., Wimpenny, J.W.T. & van Loosdrecht, M.C.M. (2001) Individual-based modelling of biofilms. *Microbiology* **147**:2897–2912. [doi:10.1099/00221287-147-11-2897](https://doi.org/10.1099/00221287-147-11-2897)
13. Lardon, L.A., Merkey, B.V., Martins, S. *et al.* (2011) iDynoMiCS: next-generation individual-based modelling of biofilms. *Environmental Microbiology* **13**:2416–2434. [doi:10.1111/j.1462-2920.2011.02414.x](https://doi.org/10.1111/j.1462-2920.2011.02414.x)
14. Luedeking, R. & Piret, E.L. (1959) A kinetic study of the lactic acid fermentation. Batch process at controlled pH. *Journal of Biochemical and Microbiological Technology and Engineering* **1**:393–412. [doi:10.1002/jbmte.390010406](https://doi.org/10.1002/jbmte.390010406)
15. Monod, J. (1949) The growth of bacterial cultures. *Annual Review of Microbiology* **3**:371–394. [doi:10.1146/annurev.mi.03.100149.002103](https://doi.org/10.1146/annurev.mi.03.100149.002103)
16. Picioreanu, C., van Loosdrecht, M.C.M. & Heijnen, J.J. (1998) Mathematical modeling of biofilm structure with a hybrid differential-discrete cellular automaton approach. *Biotechnology and Bioengineering* **58**:101–116. [doi:10.1002/(SICI)1097-0290(19980405)58:1<101::AID-BIT11>3.0.CO;2-M](https://doi.org/10.1002/(SICI)1097-0290(19980405)58:1%3C101::AID-BIT11%3E3.0.CO;2-M)
17. Pirt, S.J. (1965) The maintenance energy of bacteria in growing cultures. *Proceedings of the Royal Society B* **163**:224–231. [doi:10.1098/rspb.1965.0069](https://doi.org/10.1098/rspb.1965.0069)
18. Pirt, S.J. (1967) A kinetic study of the mode of growth of surface colonies of bacteria and fungi. *Journal of General Microbiology* **47**:181–197. [doi:10.1099/00221287-47-2-181](https://doi.org/10.1099/00221287-47-2-181)
19. Ratkowsky, D.A., Olley, J., McMeekin, T.A. & Ball, A. (1982) Relationship between temperature and growth rate of bacterial cultures. *Journal of Bacteriology* **149**:1–5. [doi:10.1128/jb.149.1.1-5.1982](https://doi.org/10.1128/jb.149.1.1-5.1982)
20. Ratkowsky, D.A., Lowry, R.K., McMeekin, T.A., Stokes, A.N. & Chandler, R.E. (1983) Model for bacterial culture growth rate throughout the entire biokinetic temperature range. *Journal of Bacteriology* **154**:1222–1226. [doi:10.1128/jb.154.3.1222-1226.1983](https://doi.org/10.1128/jb.154.3.1222-1226.1983)
21. Rittmann, B.E. & McCarty, P.L. (2001) *Environmental Biotechnology: Principles and Applications.* McGraw-Hill, New York.
22. Roels, J.A. (1983) *Energetics and Kinetics in Biotechnology.* Elsevier Biomedical Press, Amsterdam.
23. Rosso, L., Lobry, J.R. & Flandrois, J.P. (1993) An unexpected correlation between cardinal temperatures of microbial growth highlighted by a new model. *Journal of Theoretical Biology* **162**:447–463. [doi:10.1006/jtbi.1993.1099](https://doi.org/10.1006/jtbi.1993.1099)
24. Rosso, L., Lobry, J.R., Bajard, S. & Flandrois, J.P. (1995) Convenient model to describe the combined effects of temperature and pH on microbial growth. *Applied and Environmental Microbiology* **61**:610–616. [doi:10.1128/aem.61.2.610-616.1995](https://doi.org/10.1128/aem.61.2.610-616.1995)
25. Stewart, P.S. (1998) A review of experimental measurements of effective diffusive permeabilities and effective diffusion coefficients in biofilms. *Biotechnology and Bioengineering* **59**:261–272. [doi:10.1002/(SICI)1097-0290(19980805)59:3<261::AID-BIT1>3.0.CO;2-9](https://doi.org/10.1002/(SICI)1097-0290(19980805)59:3%3C261::AID-BIT1%3E3.0.CO;2-9)
26. Stewart, P.S. (2003) Diffusion in biofilms. *Journal of Bacteriology* **185**:1485–1491. [doi:10.1128/jb.185.5.1485-1491.2003](https://doi.org/10.1128/jb.185.5.1485-1491.2003)
27. Walters, M.C., Roe, F., Bugnicourt, A., Franklin, M.J. & Stewart, P.S. (2003) Contributions of antibiotic penetration, oxygen limitation, and low metabolic activity to tolerance of *Pseudomonas aeruginosa* biofilms. *Antimicrobial Agents and Chemotherapy* **47**:317–323. [doi:10.1128/aac.47.1.317-323.2003](https://doi.org/10.1128/aac.47.1.317-323.2003)
28. Wanner, O. & Gujer, W. (1986) A multispecies biofilm model. *Biotechnology and Bioengineering* **28**:314–328. [doi:10.1002/bit.260280304](https://doi.org/10.1002/bit.260280304)
29. Werner, E., Roe, F., Bugnicourt, A. *et al.* (2004) Stratified growth in *Pseudomonas aeruginosa* biofilms. *Applied and Environmental Microbiology* **70**:6188–6196. [doi:10.1128/aem.70.10.6188-6196.2004](https://doi.org/10.1128/aem.70.10.6188-6196.2004)
30. Zwietering, M.H., Jongenburger, I., Rombouts, F.M. & van 't Riet, K. (1990) Modeling of the bacterial growth curve. *Applied and Environmental Microbiology* **56**:1875–1881. [doi:10.1128/aem.56.6.1875-1881.1990](https://doi.org/10.1128/aem.56.6.1875-1881.1990)
31. Zwietering, M.H., Wijtzes, T., de Wit, J.C. & van 't Riet, K. (1992) A decision support system for prediction of the microbial spoilage in foods. *Journal of Food Protection* **55**:973–979. [doi:10.4315/0362-028X-55.12.973](https://doi.org/10.4315/0362-028X-55.12.973)
