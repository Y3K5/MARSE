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
from bulk pH by more than a unit. MARSE therefore evaluates a `ph` factor at
the pH of each voxel, which the acids, bases and ions there set (§3.8), not at
the pH of the bulk.

Its slope, which the rate Jacobian needs, is

$$
\frac{d\gamma_{\mathrm{pH}}}{d\,\mathrm{pH}} = \frac{n'd - n d'}{d^2},
\qquad
n = (\mathrm{pH} - \mathrm{pH}_{\min})(\mathrm{pH} - \mathrm{pH}_{\max}),
\quad
d = n - (\mathrm{pH} - \mathrm{pH}_{\mathrm{opt}})^2,
$$

with $n' = 2\,\mathrm{pH} - \mathrm{pH}_{\min} - \mathrm{pH}_{\max}$ and
$d' = n' - 2(\mathrm{pH} - \mathrm{pH}_{\mathrm{opt}})$, and zero outside the
limits.

Implemented as `marse.microbes.cardinal.cardinal_ph` and `cardinal_ph_slope`.

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
$\mathrm{C}_{c_j}\mathrm{H}_{h_j}\mathrm{O}_{o_j}\mathrm{N}_{n_j}$, which may
also hold phosphorus, potassium, chlorine and sodium ($p_j$, $k_j$, $l_j$,
$s_j$ atoms), with charge $z_j$, counted per mol of that formula unit (per
C-mol for biomass written per carbon atom). Its composition in the conserved
quantities is

$$
I_{j,\mathrm{C}} = c_j,
\qquad
I_{j,\mathrm{N}} = n_j,
\qquad
I_{j,e} = \gamma_j = 4c_j + h_j - 2o_j - 3n_j + 5p_j + k_j - l_j + s_j - z_j ,
$$

and $I_{j,\mathrm{P}} = p_j$, $I_{j,\mathrm{K}} = k_j$, $I_{j,\mathrm{Cl}} = l_j$,
$I_{j,\mathrm{Na}} = s_j$ for each of those elements that a network contains.

$\gamma_j$ is the **degree of reduction** (Roels 1983): the electrons released
when the component is oxidised completely to CO₂, H₂O, NH₄⁺, phosphate, K⁺,
Cl⁻, Na⁺ and H⁺, the reference state in which it is zero. Phosphorus counts
at its valence in phosphate, +5, so H₃PO₄, H₂PO₄⁻ and the salt ions all hold
no electrons, and a phosphorylated sugar holds exactly the electrons of the
sugar. Glucose has 24, biomass
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
\sum_j \nu_{pj}\, I_{jk} = 0 \qquad \text{for } k = \mathrm{C},\ \mathrm{N},\ e,
\text{ and each element present} .
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

The hydrogen balance then holds automatically. The conditions give
$\sum_j \nu_{pj}\gamma_j = 0$ with the carbon, nitrogen and element sums all
zero, so

$$
\sum_j \nu_{pj}\, h_j = \sum_j \nu_{pj}\,(2o_j + z_j) = -2\nu_{\mathrm{H_2O}} - \nu_{\mathrm{H^+}} .
$$

That is exactly the hydrogen balance with water and protons included. So every
row that satisfies the conditions is a complete chemical equation once H₂O and
H⁺ are added, and no row that fails them can be completed. This is why these
quantities suffice, and why water and protons must never be listed as
components: they hold no carbon, nitrogen, electrons or other element, so no
balance would constrain them.

**Balancing.** A process gives some coefficients and names a set $B$ of at most
as many components as there are quantities, whose coefficients $x_b$ the
balances determine:

$$
\sum_{b \in B} x_b\, I_{bk} = -\sum_{j \notin B} \nu_{pj}\, I_{jk}
\qquad \text{for every quantity } k .
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

### 3.7 Rates

A runnable network gives each process p a rate law

$$
r_p = k_p \; c_{a(p)} \prod_{i \in F(p)} f_i(c_{j(i)}),
$$

where $k_p$ is `maximum_per_h`, $a(p)$ the component the rate is proportional
to (by default the biomass of a growth process), and each factor $f_i$ is
dimensionless:

$$
\text{monod: } \frac{S}{K+S}, \qquad
\text{inhibition: } \frac{K_I}{K_I+S}, \qquad
\text{haldane: } \frac{S}{K+S+S^2/K_I}.
$$

Concentrations then change as

$$
\frac{dc}{dt} = N^{\mathsf T} r(c),
$$

with $N$ the stoichiometric matrix of §3.6. One rate drives every component a
process touches. Uptake is therefore growth divided by yield, exactly, and
every product is paid for by what is consumed (known defect 2 of the version 1
engine). All processes act on the same state at once (known defect 7).

Two rules are enforced when a network is read. First, a component limits a
process at most once. Applying a substrate's Monod term twice was known
defect 4. Second, a process may consume only components its rate depends on:
through $a(p)$, or through a monod or haldane factor. Then $r_p \to 0$ as any
reactant $c_j \to 0$, which makes the system *quasi-positive*: no
concentration can be driven below zero by the equations themselves. A
component that is genuinely never limiting may be declared
`assumed_in_excess`, the explicit form of an assumption the activated-sludge
models make for ammonium. A run in which it runs out stops with an error.

Implemented as `marse.microbes.kinetics`; the switching functions are those of
§1.3 and §1.5 (`marse.microbes.growth`).

### 3.8 Acid–base equilibria and pH

A network sets a pH when one of its components is an acid–base total
([networks.md](networks.md#acids-bases-and-ph)). Protons are never components.
In every voxel, the hydrogen ion concentration $h$ is the one that makes the
liquid electrically neutral.

**Totals.** A total $T_k$ holds every protonation state of one acid: lactic
acid and lactate, or carbonic acid, bicarbonate and carbonate. It is written in
its most protonated form, of charge $z_{k0}$, with pKa values
$pK_1 < \dots < pK_n$. These are conditional constants, for the liquid's
temperature and ionic strength. The fraction that has lost $j$ protons is

$$
\alpha_j(h) = \frac{K_1 \cdots K_j / h^j}{\sum_{i=0}^{n} K_1 \cdots K_i / h^i}
$$

(Stumm and Morgan 1996, ch. 3), so the total's mean charge is
$\bar z_k(h) = z_{k0} - \sum_j j\,\alpha_j(h)$.

**The charge balance.** With ions of fixed charge $z_s$, which are the charged
components without pKa values,

$$
F(h) = h - \frac{K_w}{h} + \sum_k T_k\, \bar z_k(h) + \sum_s z_s\, c_s = 0 .
$$

Every term rises with $h$: $d\bar z_k/dh = \sigma_k^2(h)/h$, where $\sigma_k^2$
is the variance of the number of protons lost. So $F$ rises from $-\infty$ as
$h \to 0$ to $+\infty$, and has exactly one root.

**Solving it.** Newton's method runs in $\log h$, from pH 7, inside a bracket
from pH −1 to pH 17 that every evaluation narrows. A step that would leave the
bracket, or that is not at most half the step before it, bisects the bracket
instead: the safeguarded Newton method of Press et al. (2007, §9.4). A voxel
is done when its step falls below $10^{-13}$ in $\log h$, or when $|F|$
reaches the rounding floor of its own terms, $16\,\varepsilon \sum
|\text{terms}|$.

Both safeguards were found by failures. Near the root, rounding in $F$ can be
larger than a step of $10^{-13}$ resolves, and Newton's method stepped between
two neighbouring values of $h$ for ever until the rounding floor stopped it.
Far from the root, in one plaque voxel, each Newton step landed exactly on the
other end of an unchanging bracket, between pH 3.00 and 5.76, until steps had
to halve. In plaque-like compositions the solve takes about 5 to 15
iterations.

**How h responds.** Rates that depend on pH need $\partial h/\partial c$. The
implicit function theorem gives it exactly:

$$
\frac{\partial h}{\partial c_k} = -\frac{\partial F/\partial c_k}{\partial F/\partial h},
\qquad
\frac{\partial F}{\partial h} = 1 + \frac{K_w}{h^2} + \sum_k \frac{T_k\,\sigma_k^2}{h} ,
$$

where $\partial F/\partial c_k$ is the charge a component carries at $h$: $\bar z_k$
for a total, $z_s$ for an ion. Acids and anions raise $h$, cations lower it,
and uncharged components leave it where it is. Concentrations below zero count
as zero, as they do in the rates, and have no slope.

**pH in rates.** A `ph` factor is the CPM of §2.3 at the local
$\mathrm{pH} = -\log_{10}(h / 1000\ \mathrm{mol\,m^{-3}})$. In the rate
Jacobian it responds to every charged component at once:

$$
\frac{\partial \gamma_{\mathrm{pH}}}{\partial c_k}
= \frac{d\gamma_{\mathrm{pH}}}{d\,\mathrm{pH}} \cdot \frac{-1}{h \ln 10} \cdot
\frac{\partial h}{\partial c_k} .
$$

A `dissociated` factor of a total is the number of protons each unit of it has
lost at the local $h$, $n_k(h) = \sum_j j\,\alpha_j(h)$: $K_a/(K_a + h)$ for a
single pKa. Its slope is $dn_k/dh = -\sigma_k^2/h$, and in the Jacobian it
responds to every charged component through $\partial h/\partial c_k$ in the
same way.

**Fixed charges and their counter-ions.** Groups on bacteria and in the
plaque matrix, such as carboxyl groups, are particulate totals: most of
plaque's buffering between pH 4 and 7 comes from its cell walls and matrix
(Shellis and Dibdin 1988). The cations that neutralise these groups are held
by them, not free to diffuse. They must
also be released as the groups take up protons. Held fixed, they would leave
the groups' charge unbalanced as the plaque acidified, and the lactate leaving
the plaque would carry the acid's protons with it: an oral scene built that
way turned the whole salivary film above the plaque to pH 3.3 after a sugar
rinse. Two processes keep the held cations $b$ equal to the groups' charge,
$n(h)\,T$:

$$
r_{\text{bind}} = k\,T\,n(h), \qquad r_{\text{release}} = k\,b ,
$$

with $k$ fast against every other time scale, so that $b$ stays at
$n(h)\,T$. As the groups take up the acid's protons, their cations are
released into the liquid, and lactate leaves the plaque with them, as a salt;
the acid stays on the groups until alkali from the saliva takes it off. This
is the mechanism by which the fixed buffer prolongs the low-pH phase (Dibdin
1990).
**Conservation.** Totals are components like any other, so every process still
conserves carbon, nitrogen, electrons and each element present (§3.6).
Protonation changes none of them. The protons a process releases or takes up
are whatever the charge balance then requires. They need no bookkeeping,
because $h$ is not part of the state but a function of it.

**Assumptions.**

- **Local electroneutrality.** Totals and ions diffuse independently, and $h$
  makes every voxel neutral. This treats H⁺, the fastest ion, as moving at once
  wherever neutrality needs it. It neglects the diffusion potential that
  couples the other ions. If every charged component diffuses at one
  diffusivity and the liquid is neutral by itself, as it is when fixed charges
  hold their counter-ions, diffusion carries no charge and this is exact.
  Components that diffused at different rates would separate charge, which
  would show as a pH artifact, so the oral scenes give every charged component
  one diffusivity.
- **Conditional constants.** pKa and pKw are for the liquid's temperature and
  ionic strength. Activities are not computed.
- **Instant equilibria.** Protonation is far faster than diffusion across a
  voxel or any metabolism.

Implemented as `marse.chemistry.acid_base`, and the factor in
`marse.microbes.kinetics`.

### 3.9 Oxygen roles

How a species lives with oxygen decides where in a biofilm it can grow, so it
is declared rather than left implicit in its rates
([networks.md](networks.md#oxygen-roles)). With the rate of §3.7,

$$
r_p = k_p \, c_{a(p)} \prod_i f_i(\mathbf{c}),
$$

a growth process *needs* oxygen if it consumes oxygen ($N_{p,\mathrm{O_2}} < 0$)
or has a monod or haldane factor on it, either of which is zero at
$c_{\mathrm{O_2}} = 0$. An inhibition factor $K_I/(K_I + c_{\mathrm{O_2}})$
instead equals 1 without oxygen and falls as oxygen rises. Each role constrains
the processes of its species:

| Role | Constraint |
|---|---|
| obligate aerobe | every growth process consumes oxygen |
| microaerophile | every growth process consumes oxygen, through a haldane factor $c/(K + c + c^2/K_I)$, which falls at high oxygen |
| facultative | at least one growth process consumes oxygen and at least one does not |
| aerotolerant | no process needs oxygen |
| obligate anaerobe | no process needs oxygen, and every growth process has an inhibition factor on it |

So an anaerobe's growth rate at $c_{\mathrm{O_2}} = 0$ is exactly what its
other factors give. The version 1 periodontal anaerobes, whose growth carried
a Monod term on oxygen and was therefore exactly zero without it (known
defect 3), cannot be written in version 2.

The roles state what the network models, not every capability the organism
has. They constrain the form of the rates, not their constants: whether an
inhibition constant is small enough to stop growth in air is a parameter, and
is graded like any other ([parameters.md](parameters.md)).

Implemented in `marse.schemas.network`.

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

### 4.7 Finite-volume transport on voxels

The spatial engine divides a box above a flat substratum into cubic voxels of
edge $h$, in one, two or three dimensions; the last axis is height
(`marse.spatial.grid`). Each voxel holds the mean concentration of each
component. Between two neighbouring voxels, the flux through their shared face
is

$$
F = -D\,\frac{c_{\text{upper}} - c_{\text{lower}}}{h},
$$

and a voxel changes by what enters through its faces less what leaves:

$$
\frac{dc_i}{dt} = \frac{1}{h}\sum_{\text{axes}} \left(F_{\text{lower face}} - F_{\text{upper face}}\right).
$$

Each face flux is computed once and used with opposite signs for its two
voxels. Summed over the box, every interior flux cancels, so diffusion can
change the amount inside only through the boundary. That is the property the
ledger (§9.5) checks. The boundaries are:

- **Lateral faces are periodic.** The box is one tile of a surface that
  repeats sideways.
- **The substratum is impermeable.** No flux crosses the bottom face.
- **The top face touches the bulk liquid,** held at fixed concentrations. The
  top voxel's centre is half a voxel below the face, so the flux there is
  $-D (c_{\text{bulk}} - c)/(h/2)$.

The scheme is second-order accurate in $h$. Its exact eigenvectors make that
checkable. A product of cosines that respects the boundaries (period $L$
sideways; zero slope at the substratum and zero at the top face in height)
decays at the discrete rate

$$
\lambda_h = -D \sum_{\text{axes}} \frac{4}{h^2} \sin^2\!\left(\frac{k h}{2}\right)
= -D \sum k^2 \left(1 - \frac{k^2 h^2}{12} + \dots\right),
$$

against $-D\sum k^2$ in the continuum. Around a point release, three
dimensions dilute a solute faster than two do: the Gaussian falls with
$(4\pi D t)^{3/2}$ instead of $4\pi D t$ (`point_source_diffusion_3d`). A
two-dimensional model therefore overstates how far a metabolite reaches a
neighbouring colony. The operator is second order on that release too
(docs/validation.md, "Transport in one, two and three dimensions").

A component without a diffusivity does not move: biomass stays in its voxel
until Stage 2d adds spreading. A missing lateral axis stands for a depth of
one voxel, so a 1-D column and a 3-D box report areal amounts (per m² of
substratum) in the same units.

---

### 4.8 A salivary film and the mouth

Plaque in the mouth does not sit under a well-mixed bulk liquid. It sits
under a film of saliva about 0.1 mm thick (Collins and Dawes 1987). The film
moves over the teeth at 0.8 to 7.6 mm per minute, depending on the site
(Dawes et al. 1989), and is renewed from the saliva in the mouth. The time a
substance takes to clear from plaque into the film grows with the length of
plaque the film has already crossed, and falls with the film's velocity
(Dawes 1989). Dibdin (1990) modelled the cariogenic challenge through such a
film. A domain with a `film` and a `mouth` models it the same way.

**The film** is the top $\delta$ of the box. Its surface is at the air, so
nothing crosses the top face unless the domain opens it to gases (§4.10). A
film with a free surface moves at

$$
u(\zeta) = \tfrac{3}{2}\,\bar u\,(2\zeta - \zeta^2), \qquad \zeta = \frac{z - z_0}{\delta},
$$

zero at the plaque and fastest at its surface. Its shear at the plaque is
$3\bar u/\delta$, which §6.4 uses for the delivery of cells. At a site the
film reaches after crossing a length $l$ of plaque, each layer is replaced by
fresh saliva from the mouth at

$$
k(\zeta) = \frac{u(\zeta)}{l},
\qquad
\frac{\partial c}{\partial t} = \dots + k(\zeta)\,(m - c)
$$

for every dissolved component, where $m$ is the pool's concentration. This
is a well-mixed film renewed at the rate of plug flow. For plaque upstream
that behaves like the site, it has the same steady state as the film leaving
the plaque.

**Over a plaque that spreads** (§6.1), the plaque's surface $L$ moves: it
rises as the plaque grows and falls when plaque is brushed off. The film
rides on it. The share of each voxel above $L$ is liquid, and it is renewed
at $u(\zeta)/l$ at the liquid's centre, with $\zeta = (z - L)/\delta$ at most
one. Liquid more than a film's thickness above the surface, which the box
leaves above a thinned plaque, moves with the film's surface; it stands for
saliva the film's surface renews, and its volume, at most the plaque brushed
off, is small next to the mouth's (63 µm against the 3,750 µm the pool holds
over the plaque at rest in the scenes). The renewal is recomputed at the
start of every span from the plaque's height then. On a plaque as high as the
film's underside it is the renewal above, exactly, and food left on the teeth
(§4.9) is placed in the film's thickness above the surface.

**The mouth** follows Dawes's (1983) model of sugar clearance:

- **Filling.** The liquid in the mouth grows from its resting volume, RESID,
  at the salivary flow $Q$, until it reaches VMAX (Lagerlöf and Dawes 1984
  measured both).
- **Swallowing.** A swallow is an incomplete syphon. It takes the volume back
  to RESID and leaves every concentration as it was.
- **Stimulation.** The flow is $Q = Q_u + Q_s\, s/(K + s)$, where $s$ is the
  concentration of a stimulus, such as sugar, in the mouth.
- **Secretion.** The glands secrete resting saliva at $Q_u$ and saliva closer
  to stimulated saliva as the flow rises. Bicarbonate, for one, rises steeply
  with flow (Bardow et al. 2000).

RESID and VMAX count the liquid on every surface of the mouth, and the film
over the modelled plaque is part of it. That film is in the box. So the pool
the box exchanges with has volume $V - A\delta$, where $A$ is the plaque area
the box stands for, and a thickness over it of $H = (V - A\delta)/A$. Per unit
area of plaque, the pool's amount $n = Hm$ of each dissolved component
changes as

$$
\frac{dn}{dt} = c_{\text{secreted}}(Q)\,\frac{dH}{dt} - \sum_{v} k_v\,(m - c_v)\,a ,
$$

where the sum runs over the voxels of the box and $a$ is each voxel's volume
per unit area of substratum. At a swallow, $H$ falls from $H_{\max}$ to
$H_{\text{resid}}$ and $n$ falls in proportion. Between swallows, the pool
dilutes as it grows. After $N$ swallows at constant flow, with nothing taken
up, a component has fallen to $(H_{\text{resid}}/H_{\max})^N$ of what it was.

**Assumptions.**

- **The pool is well mixed**, and the box stands for all of the plaque area
  $A$.
- **Each layer of the film is renewed independently**, at $u/l$. Within a
  voxel the film is mixed.
- **The flow over each span** (§9.9) follows the stimulus the pool holds at
  the span's start, diluted by secretion and added to by what an intake
  brings (§4.9) and by what the plaque gave back over the span before. A
  change in what the plaque takes up or gives back reaches the flow one span
  later, at most a minute. This changes when the swallows fall, not what is
  conserved.

Implemented as `marse.oral` and `marse.core.reservoir`.

### 4.9 The diet

After sugar, the pH of plaque falls within minutes to a minimum, then returns
over half an hour or more: the Stephan curve (Stephan 1944). How deep and how
long depends on how the sugar comes. A high-sugar eater differs from a
low-sugar eater in more than the amount: sugar stays in the mouth for longer,
sipped in drinks or sucked from sweets, and food retained on the teeth holds
it where the plaque is. Particles of starchy foods stayed on the teeth for up
to 20 minutes and accumulated sugars and acids (Kashket, Zhang and Van Houte
1996). A `diet` lists intakes, each from a start for a duration $T$, one at a
time. In the notation of §4.8:

- **A rinse** of volume $V_r$ and composition $c_r$ is taken in at once: $V$
  gains $V_r$, and $n$ gains $c_r V_r / A$. The mouth holds it without
  swallowing, secreting as it tastes it, and at the end expels everything
  above RESID. $H$ falls to $H_{\text{resid}}$ and $n$ in proportion, as at a
  swallow.
- **A drink** of volume $V_d$ and composition $c_d$ flows in at $q = V_d/T$:
  $dV/dt = Q + q$, and $dn/dt$ gains $c_d\, q/A$. The mouth swallows at VMAX
  as before.
- **A food** releases amounts $R$ at $R/T$ without liquid: $dn/dt$ gains
  $R/(AT)$.
- **Chewing.** A food that is chewed, as gum or a meal is, adds a flow of its
  own while it lasts: $Q = Q_u + Q_c + Q_s\, s/(K + s)$. Dawes and
  Macpherson (1992) found that flavoured gum raised the flow to 10 to 12
  times the resting 0.47 mL per minute in the first minute, falling within
  about 10 minutes to the flow gum base alone gave, and Dawes and Kubieniec
  (2004) found it still above the resting flow after two hours. $Q_c$ is
  that lasting part. The glands secrete saliva closer to stimulated saliva
  as the flow rises, whatever raises it. In the model of Dibdin, Dawes and
  Macpherson (1995), sugar-free gum chewed during a cariogenic challenge
  raised the pH of plaque quickly.
- **Mixing.** While an intake lasts, each voxel of the film also exchanges
  with the pool at $k_{\text{mix}}$, so that $k(\zeta) = u(\zeta)/l +
  k_{\text{mix}}$. Dibdin (1990) mixed the film with the mouth's liquid while
  sugar was in the mouth. The default, once a second, brings the film to the
  mouth's composition within seconds.
- **Food left on the teeth.** An amount $M$ per unit area of a particulate
  component is placed in the film over its region when the intake ends, at
  $M/\delta$ in each voxel. It does not move and is not renewed. A process of
  the network, such as $r = k\,c_{\text{food}}$, releases what dissolves from
  it, which diffuses into the plaque and out through the film.

What an intake brings is booked as eaten on the ledger of the box and the
mouth together (§9.5), and what a rinse takes away as expelled. Food placed in
the film enters the box, and the box's own ledger books it so.

**Assumptions.**

- **Steady intakes.** A drink or a food enters at a constant rate over its
  duration. Sips and bites are smoothed out; a series of short intakes can
  resolve them.
- **A food brings no liquid.** A dissolving sweet adds little volume beside
  the saliva it stimulates.
- **Chewing adds a constant flow.** The first minutes of flavour, when the
  flow is several times higher, are left out; a short intake that releases a
  stimulus can stand for them.
- **Mixing only while an intake lasts**, the same through the film's depth.
  Between intakes the film is renewed by its flow alone.
- **Retained food lies in the film**, above the plaque, not within it.

**What the oral scenes reproduce.** With a plaque buffer whose groups release
their cations as they take up protons (§3.8), the scenes of
`examples/environments/oral`, calibrated to a Stephan curve, also move the way
the literature reports when one thing changes. Sipping and food left on the
teeth keep the plaque below pH 5.5 for longer. Low salivary flow gives a
greater fall that stays low for longer (Lingström and Birkhed 1993). A slower
film gives a lower minimum and a slower return (Macpherson and Dawes 1991).
Thicker plaque stays acid for longer, although Dawes and Dibdin (1986) found
an optimum thickness for the lowest pH. More fixed buffer makes the fall
shallower and the low phase longer (Dibdin 1990). The numbers are in
[validation.md](validation.md#the-stephan-curve).

Implemented as `marse.oral.diet`, with the measures of a curve in
`marse.oral.stephan`.

### 4.10 Oxygen from the air

Plaque takes its oxygen from the air through the salivary film, and how far
the oxygen reaches sets where aerobes, facultative anaerobes and anaerobes
can live. Under saliva, plaque became anoxic below about 220 µm, and below
about 150 µm after a rinse of sucrose (von Ohle et al. 2010). A domain with
`air` opens the top face of the box to the gases it lists, and to nothing
else.

**At the face**, the air holds each gas at its saturation $C_s$, its
concentration in water in equilibrium with the air: 0.210 mol m⁻³ for oxygen
at 37 °C under air at one atmosphere (§4.4). In the finite-volume scheme of
§4.7, a face held at a concentration exchanges with the voxel below it across
half a voxel of height $h$:

$$
J = \frac{2D}{h}\,(C_s - c_{\text{top}})
$$

per unit area. MARSE runs it as two processes in the top layer, the air
giving at $k C_s$ and taking back at $k c_{\text{top}}$, with $k = 2D/h^2$.
The limiter keeps the voxel from going negative, as it does for any process,
and the ledger counts what crosses as an import, exchanged with the air. A
slab filling from its surface (Crank 1975), and the profile under zero-order
uptake $k_0$,

$$
c = C_s\left(1 - \frac{x}{\delta}\right)^{2}
\quad\text{within}\quad
\delta = \sqrt{\frac{2 D C_s}{k_0}}
$$

of the surface, are both followed at second order in the voxel
([validation](validation.md#oxygen-from-the-air)).

**Under a film**, the face is the film's surface. The mouth's own liquid is at
the air too. It lies over the mouth's surfaces as a film about as thick,
which comes to equilibrium with the air at the rate of its slowest mode,
$\pi^2 D / 4\delta^2$: within 1.5 s for oxygen in a film 0.1 mm deep. A
stated rate would only be a stiff number. So the pool holds each gas at $C_s$,
and the film is renewed with saturated saliva, $m = C_s$ in §4.8. Holding it
takes or gives $C_s\,dH - dn$, the difference between what the pool would
hold and what it does hold, and that too is booked as exchanged with the
air, on the ledger of the box and the mouth.

**Without a film**, the box borders the air instead of a bulk liquid, as a
colony biofilm grown on a membrane does.

**Assumptions.**

- **The air holds each gas at a stated saturation**, the same over the film
  and in the mouth. Breath changes the air in the mouth: exhaled air holds
  less oxygen than the room's.
- **Carbon dioxide stays closed.** It belongs to the carbonate total of §3.8,
  whose CO₂ the buffer calibrated in Stage S1 keeps; a gas must be a neutral
  component and not an acid-base total. Saliva that loses CO₂ to the air
  rises in pH, which is why Bardow et al. (2000) collected it under oil;
  MARSE leaves that loss out.
- **The air offers no resistance**: the face is at saturation, however fast
  the box takes the gas up.

Implemented as `marse.spatial.air` and, under a film, `marse.core.reservoir`.

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

The approximation fails during fast transients, which are precisely the
perturbation experiments MARSE is built for, so it cannot be silently trusted.
MARSE's spatial engine for configuration schema version 2 therefore does not
make it. Solutes and biomass are integrated together in time by an L-stable
implicit method (§9.8):

- with long steps it lands on the quasi-steady profile;
- with short steps it follows a transient;
- no switch between the two regimes is needed.

The ratio $\tau_D/\tau_{\text{growth}}$ is still computed, from the box height
and the fastest diffusivity and rate, and recorded in every manifest
(`time_scales`). It tells a reader how stiff the run was, not whether its
results can be trusted.

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
individual-based mechanical relaxation (Kreft et al.). They mix lineages to
different degrees, which changes conclusions about competition
([modeling-landscape.md §2](modeling-landscape.md#2-the-biomass-spreading-decision)),
so MARSE makes the mechanism a declared provider, starting with two on its
voxels ([the Stage 2d plan](stage-2d-plan.md)).

**The constraint.** Each particulate component $j$ that takes up room has a
packing density $\rho_j$, its concentration when it alone fills space. The
fraction of a voxel its biomass fills is

$$
\phi = \sum_j \frac{c_j}{\rho_j} \le 1 ,
$$

every species counted together, so the room is shared by construction.

**The displacement.** Where the biofilm is full, growth makes volume at the
rate $s = \sum_j r_j / \rho_j$, with $r_j$ the net production of component $j$
by the reactions, and the biomass flows away from it as an incompressible
material through a porous medium does, $\mathbf{u} = -\nabla p$ with
$\nabla\cdot\mathbf{u} = s$ (Alpkvist and Klapper 2007). All the components at a
point move with one velocity. In a column full from the substratum up, this
is Wanner and Gujer's displacement velocity,

$$
u(z) = \int_0^z \sum_j \frac{r_j}{\rho_j}\, dz' .
$$

With uniform growth at rate $\mu$, $u = \mu z$, so every point of the film
moves from $z_0$ to $z_0 e^{\mu t}$. Fed from the liquid instead, only a layer
as deep as the substrate penetrates grows, so a deep film thickens at the
constant rate $Y J / \rho$, with $J$ the substrate flux into it: the
one-dimensional form of the linear expansion of §6.2.

**Two mechanisms on voxels.** Both keep $\phi \le 1$, carry each component
that moves in the proportions it is found in, and conserve each exactly. Both
run in a column for now.

- **Packed** (`displacement_1d_v1`, `marse.biofilm.spreading`): the plaque of
  the oral scenes. After every step of the integration, the solid of voxel $k$
  occupies the cumulative volume between $V_{k-1}$ and
  $V_k = \sum_{i \le k} \phi_i$, in voxels, and new voxel $j$ takes whatever
  lies between $j$ and $j + 1$, each component in the proportions of the voxel
  it came from. Each component's amount below a given volume is a
  non-decreasing, piecewise-linear function of that volume, and the new voxels
  hold its differences, so nothing goes negative and everything is accounted
  for. The front is sharp: full voxels, then at most one partly filled. Where
  growth slows below a voxel's capacity, the column settles.
- **Continuum** (`continuum_pressure_v1`, `marse.spatial.spreading`): growth
  runs in place for an interval, and the excess is then pushed on by the
  pressure it makes, keeping material in order. A voxel that empties stays
  empty, and the method extends to two and three dimensions (§9.10).

What the packed plaque adds:

- **Detachment at a maximum height.** Solid pushed above the plaque's
  maximum height $L_{\max}$ leaves the column: in the mouth, the film's flow
  carries it off, and the mouth swallows it.
- **Accuracy.** Solid grows for a step before it is packed, so what passes
  $L_{\max}$ is cut a step late, a first-order error that shrinks with the
  step. A film of $L_0 = 40$ µm growing at $\mu = 0.2$ h⁻¹ under a maximum of
  60 µm detached $1.0 \times 10^{-3}$ more than the continuous solution over
  six hours.
- **The closed form.** A film growing at a constant specific rate $\mu$ and
  worn at a constant velocity $u_w$ (§6.3) thickens as
  $dL/dt = \mu L - u_w$, so
  $L(t) = u_w/\mu + (L_0 - u_w/\mu)\,e^{\mu t}$. The fixed point
  $u_w/\mu$ is unstable: below it the film wears away, above it the film
  grows. MARSE follows this as closely as the integrator's tolerance asks:
  within $9.7 \times 10^{-7}$ at a relative tolerance of $10^{-6}$
  ([validation](validation.md#plaque-that-spreads)).

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

**Wear.** In the mouth, the tongue and cheeks rub the surface of plaque.
MARSE's first closure is the simplest: a constant wear velocity $u_w$, which
takes $u_w\,\Delta t$ of solid off the top of the column over a step. Half
comes off before the step and half after (Strang splitting), which keeps the
step second order. The integrator's error control cannot see a splitting
error, so this matters. A film of 20 µm against the closed form of §6.1:

| Relative tolerance | Half the wear before the step, half after | All of it after |
|---|---|---|
| $10^{-5}$ | $9.7 \times 10^{-6}$ | $1.1 \times 10^{-2}$ |
| $10^{-6}$ | $9.7 \times 10^{-7}$ | $3.6 \times 10^{-3}$ |
| $10^{-7}$ | $9.9 \times 10^{-8}$ | $1.2 \times 10^{-3}$ |

With growth limited by what diffuses in, the surface grows more slowly as
plaque thickens, so a constant wear velocity gives a stable thickness. With
unlimited growth it cannot (§6.1).

**Brushing and flossing** take a share $f$ of the plaque off from its surface
down at a stated time. Over 59 papers and 212 brushing exercises, a brushing
with a manual toothbrush removed 42% of plaque. Depending on the plaque index
used, the figure runs from 30% to 53% (Slot et al. 2012). MARSE takes 0.42
for a brushing unless told otherwise, and asks a flossing for its share.
The same share of any food left on the teeth goes with it. In the mouth, what
a cleaning takes is expelled, and both balances count it.

### 6.4 Attachment to surfaces

A biofilm starts when cells reach a surface from the liquid, bind to it, and
stay. MARSE follows each attaching species through these steps at every face
of the substratum, the bottom face of the box (§4.7). The scenes built on this
are described in [`docs/environments.md`](environments.md).

**Delivery.** A cell of diameter $d$ is a Brownian sphere, with the
Stokes–Einstein diffusivity

$$
D = \frac{k_B T}{3\pi\eta d}.
$$

A streptococcus 0.9 µm across diffuses at 0.73 µm² s⁻¹ in water at 37 °C,
about a thousandth of the rate of a small solute. Close to a wall, any flow is
a simple shear, $u = \dot\gamma y$. A wall that captures the cells reaching it
depletes the liquid above it, in a layer that thickens downstream. Lévêque
(1928) solved this problem at steady state. At a distance $x$ downstream of
where capture begins, the concentration depends on the height $y$ and on $x$
only through $\eta = y\,(\dot\gamma/9Dx)^{1/3}$, as
$c = c_\infty \int_0^\eta e^{-s^3}\,ds \,/\, \Gamma(4/3)$. The flux onto the
wall is then

$$
j = k\,c_\infty, \qquad
k = \frac{1}{\Gamma(4/3)\,9^{1/3}} \left(\frac{D^{2}\dot\gamma}{x}\right)^{1/3}
= 0.5384 \left(\frac{D^{2}\dot\gamma}{x}\right)^{1/3}.
$$

In colloid deposition this is the Smoluchowski–Levich approximation, and
flow-chamber experiments on bacteria are analysed with it (Busscher and van der
Mei 2006). The transfer velocity $k$, multiplied by the cells per volume in the
bulk, gives the cells arriving per area and time.

A parallel-plate flow chamber sets $\dot\gamma$ directly. A film of liquid
flowing over a surface with its top free, such as saliva on a tooth, has a
half-parabolic profile, $u = \tfrac32 \bar u\,(2y/\delta - y^2/\delta^2)$ for a
film of thickness $\delta$ and mean velocity $\bar u$. Its wall shear rate is
therefore $3\bar u/\delta$.

**Binding, limited by the free area.** A fraction $\alpha$ of the delivered
cells binds. This attachment efficiency belongs to the species on that
material, under its conditioning film. Meinders, van der Mei and Busscher
(1995) define a deposition efficiency the same way: the measured initial
deposition rate over the rate the Smoluchowski–Levich flux predicts.

Bound cells block the area around them. A species $s$ bound at $n_s$ cells per
area, each blocking an area $a_s$, reduces binding in proportion to the blocked
fraction:

$$
J_s = \alpha_s\,k_s\,c_s\,B(\theta), \qquad
B(\theta) = \max\!\left(0,\; 1 - \frac{\theta}{\theta_J}\right), \qquad
\theta = \sum_s n_s a_s .
$$

$\theta_J = 0.547$ is the jamming limit of random sequential adsorption (Feder
1980). Disks placed at random, each kept only if it overlaps none already
there, stop fitting once they cover 54.7% of a plane. With $B$ linear in
$\theta$, this is the blocked-area model that flow-chamber deposition is
fitted with (Busscher and van der Mei 2006). A cell blocks more than its own
footprint, because others cannot bind in a neighbourhood around it. Every
attaching species on a face competes for the same free area.

**Reversible, then locked.** A cell binds reversibly at first, into the
species' reversible component $R$. From there it detaches back into the liquid
at the rate $k_{\mathrm{off}}$, or locks into the species' biomass $L$ at the
rate $k_{\mathrm{lock}}$:

$$
\frac{dR}{dt} = J - (k_{\mathrm{off}} + k_{\mathrm{lock}})\,R, \qquad
\frac{dL}{dt} = k_{\mathrm{lock}}\,R + \text{growth}.
$$

Locking stands for the strengthening of adhesion over the first minutes of
contact. On saliva-coated enamel, the adhesion force of streptococci rose from
−0.7 to −10.3 nN as contact lasted from 0 to 120 s (Mei et al. 2009). Locked
cells are biomass and grow through the network's growth processes.
Reversibly bound cells last minutes and take part in no process, which the
schema checks.

**The closed form.** One species with no growth and $R + L < n_J = \theta_J/a$
follows a linear system, with $j_0 = \alpha k c$:

$$
\frac{d}{dt}\begin{pmatrix} R \\ L \end{pmatrix}
= \begin{pmatrix} -(j_0/n_J + k_{\mathrm{off}} + k_{\mathrm{lock}}) & -j_0/n_J \\
k_{\mathrm{lock}} & 0 \end{pmatrix}
\begin{pmatrix} R \\ L \end{pmatrix}
+ \begin{pmatrix} j_0 \\ 0 \end{pmatrix}.
$$

On the boundary $R + L = n_J$, the total can only fall,
$\tfrac{d}{dt}(R + L) = -k_{\mathrm{off}} R \le 0$, so the solution never
leaves the region where the system is linear. With locking, every cell ends up
locked and the surface jammed, $(R, L) \to (0, n_J)$, along the two decaying
modes of the matrix. Without locking, $R$ settles on the Langmuir balance of
binding and detachment,

$$
n_\infty = \frac{j_0}{j_0/n_J + k_{\mathrm{off}}}.
$$

Implemented as `marse.validation.analytical.adhesion_kinetics`: the reference
the engine is tested against.

**Conservation at the substratum.** In the engine, bound cells live in the
bottom layer of voxels. An areal density is a concentration times the voxel
height $h$, and a number of cells is an amount divided by the carbon in one
cell, $q$. Deposition enters each bottom voxel through its substratum face,
and detachment leaves through it, as a face flux in the sense of §9.8:

$$
F_{\mathrm{sub}} = \alpha\,k\,c\,q\,B(\theta) - k_{\mathrm{off}}\,R\,h .
$$

- **Counted, not inferred.** The ledger books this flux as an import, beside
  the exchange through the top face.
- **Locking conserves exactly.** It converts $R$ into $L$, two components the
  schema requires to have the same formula. So it conserves carbon, nitrogen
  and electrons exactly. The engine runs it as an extra row of the
  stoichiometric matrix, one per species, acting in the bottom layer only.
- **An analytic Jacobian.** Blocking couples every attaching species on a
  face. Like the rate laws, the exchange is evaluated at concentrations
  clamped at zero and differentiated consistently, so nothing responds below
  zero. Past jamming, binding has no slope.
- **Limited like any transfer.** The limiter of §9.8 counts deposition as an
  arrival and detachment as a departure, which it scales like any other.

**What this model leaves out.**

- **Sedimentation.** Cells denser than the liquid also settle. The scenes this
  serves have surfaces that are vertical or face down, so it is omitted.
- **Charge and surface energy.** $\alpha$ is stated per species and material.
  The next increment of the environments programme computes it from the
  extended DLVO theory, for surfaces and waters without adhesion data.
- **Coadhesion.** Cells bind to the substratum, not to cells already bound.
  Once growth covers the surface past jamming, binding stops.
- **A fixed suspension.** Binding does not deplete the cells in the bulk
  liquid, and cells are not tracked through the liquid inside the box.
- **Roughness.** Surfaces rougher than about $R_a = 0.2$ µm retain more
  plaque (Bollen et al. 1997), mainly by sheltering cells from removal. That
  belongs to retention, not to binding.

Implemented across four modules:

- `marse.spatial.colloids`: delivery;
- `marse.spatial.surface`: which material each face is;
- `marse.microbes.adhesion`: binding, detachment and locking;
- `marse.core.implicit`: the exchange within the integration.

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

Milliseconds, for a process to be simulated over days. In three dimensions on
2 µm voxels the limit is a quarter of a millisecond. Explicit stepping of the
full coupled system is therefore not viable. The version 1 engines avoid it by
solving the solute field to steady state (§9.2). The version 2 spatial engine
instead integrates transport implicitly (§9.8), and `marse check` prints the
explicit limit it avoids.

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

For reaction networks the check is a ledger over the three conserved
quantities of §3.6. After every step, for each quantity k,

$$
\text{residual}_k = M_k(t) - M_k(0) - \text{imports}_k + \text{exports}_k,
\qquad M_k = \sum_j I_{jk} c_j .
$$

A closed box has no imports or exports. In space, what crosses the top face
is booked from the face transfers the integrator applied (§9.8), independently
of the amounts in the box, so the check compares two separately computed
quantities. A residual beyond $10^{-9}$ of the total content
$\sum_j |I_{jk}| c_j(0)$, plus everything that has crossed, stops the run
with a `ConservationError`. That threshold sits far above rounding, which is about
$10^{-15}$ per step, and far below any real leak. The largest residual of every
run is recorded in its manifest. Implemented as `marse.core.ledger`.

### 9.6 Positive, conservative integration

Reaction networks are integrated with Heun's method, written as the average of
a state and two forward-Euler steps (Shu and Osher 1988):

$$
y^{(1)} = E_h(y^n), \qquad y^{n+1} = \tfrac12\left(y^n + E_h(y^{(1)})\right).
$$

Each Euler step is *limited process by process*. With extents
$\xi_p = h\,r_p(y)$, a species $j$ is overdrawn when
$\sum_p \max(-N_{pj},0)\,\xi_p > y_j$. Every process consuming an overdrawn
species is scaled by that species' share $y_j / \text{consumption}_j$, and a
process is limited by the scarcest species it consumes. The consequences:

- **Conservation.** Each process row conserves carbon, nitrogen and electrons
  (§3.6), and a scaled row still does. The update therefore keeps every
  balance to rounding, whatever $h$.
- **Positivity.** After limiting, no species loses more than it held, so
  every Euler step leaves it non-negative, and $y^{n+1}$ is an average of
  non-negative states: the strong-stability-preserving property of this form
  of Heun's method (Gottlieb, Shu and Tadmor 2001). Nothing is clipped.
- **Selectivity.** Only the processes that consume a depleted species slow
  down.

  A single factor for the whole increment, as in the BBKS schemes (Bruggeman
  et al. 2007), is also positive and conservative. But then one depleted
  species slows every process, and a species that starts at zero can halt the
  run. Both effects were measured in MARSE's own prototypes before this
  scheme was chosen.

Production within a step is not counted toward what may be consumed.
Counting it would let a cycle of processes pass arbitrarily large amounts
through a species within one step, and the cancelling flows would cost
conservation its precision. That was also measured, as a drift of $10^{-8}$.

Without limiting the method is second order, and the test suite measures the
order. The Euler result and $y^{n+1}$ form an embedded pair. The substep is
adapted so that, for every species,

$$
|y^{n+1}_j - y^{(1)}_j| \le \text{atol} + \text{rtol}\,\max\left(|y^n_j|, |y^{n+1}_j|, \hat y_j\right),
$$

where $\hat y_j$ is the largest value species $j$ has reached. A species that
has run down to a trace is thus held to the accuracy of its own scale, not
resolved ever more finely once it no longer matters. The substeps are chosen
deterministically, so a run replays bit for bit. Implemented as
`marse.core.integrators`.

### 9.7 Multigrid for implicit systems

An implicit step on a grid (§9.1 explains why transport must be implicit)
needs, at every Newton iteration, the solution of

$$
x - a\left(L x + B x\right) = b,
$$

where

- $L$ is the diffusion operator of §4.7 with the bulk held at zero, as a
  correction equation requires;
- $B$ holds one small matrix per voxel, the Jacobian of the reactions there,
  coupling the components within a voxel;
- $a$ is the step weight.

With the steps an implicit scheme exists to take, $aD/h^2$ reaches $10^4$ or
more. Simple iterations then remove the rough part of the error quickly but
the smooth part only very slowly.

Multigrid (Briggs, Henson and McCormick 2000) treats each part on the grid
where it is rough:

- **Smoothing.** A few sweeps of red–black Gauss–Seidel damp the rough error.
  Each voxel's components are solved together through their small block, so
  strong local chemistry does not slow the sweep.
- **Coarse-grid correction.**
  - The remaining residual is averaged onto a grid with half as many voxels
    per axis, where the smooth error is rough again.
  - It is corrected there, recursively.
  - The correction is interpolated back linearly between voxel centres.
  - Coarse operators are re-discretised, and the reaction blocks are averaged
    over each voxel's children.
- **The coarsest grid** is solved exactly, through its inverse, which MARSE
  computes by a blocked LU factorisation with partial pivoting. LAPACK would
  be faster, but OpenBLAS factorises a matrix of more than about 100 unknowns
  in parallel, and its last bits then depend on the number of threads. The
  blocked factorisation leaves those bits to matrix products, which OpenBLAS
  computes the same way on any number of threads, and to numpy's elementwise
  arithmetic.

Each cycle then reduces the error by a similar factor on any grid, measured at
0.06–0.09 from 8³ to 64×64×32 voxels. The work per cycle is proportional to
the number of voxels. One cycle preconditions a restarted GMRES (Saad and
Schultz 1986), which keeps convergence robust when the reaction blocks make
the system non-symmetric. Residuals are weighed on the same per-entry scale
the error control uses, because biomass near $10^3$ mol m⁻³ and oxygen near
$10^{-1}$ must be solved to the same relative accuracy. Every operation runs
in a fixed order, so the same system gives the same answer, bit for bit, on
any number of threads.
Implemented as `marse.spatial.multigrid`.

**A column is solved directly.** In one dimension the system couples each
voxel only to its two neighbours. Block Gaussian elimination from the
substratum up, then back substitution (the Thomas algorithm, in blocks),
solves it exactly in $O(nJ^3)$. It is three times faster than multigrid on the
1-D example, and its results agree with multigrid's to about $10^{-14}$.
Implemented as `marse.spatial.column`.

### 9.8 Implicit reaction–transport integration

In space, every voxel changes by diffusion (§4.7) and by reaction (§3.7) at
once:

$$
\frac{dc}{dt} = \nabla_h\!\cdot F(c) + N^{\mathsf T} r(c).
$$

Three ways of stepping this fail:

- **Explicit steps** must stay below the limit of §9.1, a quarter of a
  millisecond on 2 µm voxels.
- **Splitting reaction from transport** at practical steps starves the
  biofilm. Each step refills the box once, where a real biofilm is fed by a
  continuous flux.
- **A quasi-steady solute field** (§5.3) is wrong during transients, which
  are what MARSE is built to study.

MARSE therefore integrates everything together, with the two-stage, L-stable,
stiffly accurate diagonally implicit Runge–Kutta method of Alexander (1977),
$\gamma = 1 - 1/\sqrt 2$:

$$
Y_1 = y^n + \gamma h\, f(Y_1), \qquad
Y_2 = y^n + h\left((1-\gamma) f(Y_1) + \gamma f(Y_2)\right).
$$

Its stability function $R(z) = (1 + (1-2\gamma)z)/(1-\gamma z)^2$ tends to zero
as $z \to -\infty$. Stiff modes are therefore damped within one long step, and
a long step lands on the quasi-steady profile, while a short step follows a
transient at second order.

**Conservation for any solver tolerance.** The new state is not taken from
$Y_2$. It is rebuilt from what the stages imply: the face transfers
$T = h\sum_i b_i F(Y_i)$ and the process extents
$\xi = h\sum_i b_i r(Y_i)$,

$$
y^{n+1} = y^n + \nabla_h\!\cdot T + N^{\mathsf T}\xi .
$$

Each transfer leaves one voxel and enters its neighbour, and each process row
conserves carbon, nitrogen and electrons (§3.6). So the balance holds to
rounding however loosely the linear systems are solved. The transfers through
the top face are what the ledger (§9.5) books as the box's exchange with the
bulk liquid.

**Positivity.** Positivity for every step size and second order cannot be had
together: a scheme that keeps the heat equation positive at every step size is
at most first order (Bolley and Crouzeix 1978). Here $R(z)$ turns negative below
$z = -1/(1-2\gamma) \approx -2.4$, so a nearly exhausted species can undershoot.
Where $y^{n+1}$ would be negative, the limiter scales down what leaves that
voxel: its outgoing transfers and the processes consuming there, each as a
whole. It counts what arrives within the same step, and repeats until nothing
is negative. The stock-only limiter of §9.6 is the guaranteed fallback. Nothing
is clipped, and every balance survives the scaling. The update itself is
summed from what arrives in each voxel and what leaves it, kept apart, so a
voxel that loses nothing can only gain, even in rounding. Summed as a
divergence instead, it left traces of lactate of $10^{-323}$ mol m⁻³ negative
on the first steps of the 3-D example, with nothing leaving them to scale.
Below the smallest normal number, $2.2\times10^{-308}$, no relative margin
survives rounding. Scaled by 4/6, a transfer of one unit in the last place
rounds back up to one, so a voxel holding four such units and sending one
through each of its six faces cannot be scaled into balance. On a 64 × 64 × 32
box the limiter stalled on such traces and its fallback failed. A voxel holding
less than the smallest normal number therefore stops giving, instead of being
scaled. A last pass after the fallback does the same for any voxel that rounding
still leaves below zero. A voxel that gives nothing can only gain, so the
guarantee holds in floating point, not only in exact arithmetic.

**Newton's method.** Both stages are solved by Newton's method with the
analytic Jacobian of the rates (`rate_jacobian`). The rate laws have a kink at
zero, where a Monod term switches off. So Newton is projected:

- a value well above the absolute tolerance may fall at most to a tenth of
  itself per iteration, which keeps Newton on the smooth side while a large
  transient is resolved;
- a value already below the tolerance may go negative, because the stage
  solution itself can undershoot there.

Two alternatives were built and measured first:

- damping every value to stay positive fails when the stage solution is
  truly negative;
- a line search alone crawls through an oxygen front, one voxel per
  iteration.

The Jacobian is that of the rates as they are evaluated, at concentrations
clamped at zero, so a negative value has no slope. Given the slope at zero
instead, Newton's method is told that consumption still responds below zero.
On a ten-minute step from air-saturated biomass it then crept towards a
negative stage value by 4% per iteration and had not converged after 40
iterations. With the true slope it converged in 18.

A Newton solve that has not converged after 20 iterations rejects the step,
which is retried at a quarter of its length. One matrix, with the Jacobian at
the start of the step, serves both stages, because they share $\gamma h$. It is
rebuilt only when Newton converges slowly, and its systems are solved by
multigrid (§9.7).

**Error control.** The first-order result $y^n + h f(Y_1)$ differs from
$y^{n+1}$ by $h\gamma\,(f(Y_2) - f(Y_1))$. That difference is filtered through
$(I - \gamma h J)^{-1}$, so that the stiff components the method damps do not
dominate it (Hairer and Wanner 1996, §IV.8), and the filtered difference sets
the step. Errors are measured against each component's largest value
anywhere in the box. A trace of lactate seeping into the liquid far from the
colonies is then held to the accuracy lactate needs where it matters.
Measured per voxel, the first three minutes of the column below took 1,292
steps; measured per component, the first fifteen minutes take 388, with an
error of $1.0\times10^{-6}$ of each component's peak. What remains of the
start-up is real: lactate appearing from nothing, then oxygen running out at
the base within seconds, a transient spread over five decades of time. A
second-order method needs about 60 steps per decade to resolve it at a
relative tolerance of $10^{-4}$, and about 20 at $10^{-3}$. Filtering the
estimate a second time was measured too. It saved 1% of the steps and cost
50% more time, so MARSE filters once.

**The first step.** Colonies placed in fresh liquid start far from their
quasi-steady state, and a first step as long as the recording interval only
fails. The first step of a run therefore comes from the starting-step
algorithm of Hairer, Nørsett and Wanner (1993, §II.4), for order 2 and in the
error norm above. It uses the rates and their change over a trial Euler step,
so that the increment is a hundredth of the state and the local error a
hundredth of the tolerance. On the column below it proposes 0.1 ms, and the
first fifteen minutes see two rejected steps. Starting with the
fifteen-minute recording interval instead cost ten, four of them failed
Newton solves.

**Measured.** On a column of 100 voxels of 2 µm, with the example chemistry
and 800 and 400 C-mol m⁻³ of heterotroph and fermenter in the lower 100 µm:

- carbon, nitrogen and electrons are conserved to $2.4\times10^{-15}$ of their
  totals over $10^4$ steps, counting imports;
- the observed order rises from 1.74 to 1.92 as the step falls from 90 s to
  2.8 s. Below 2 at long steps is the order reduction expected of a method of
  stage order 1 on a stiff problem (Prothero and Robinson 1974);
- the actual error is 20 to 50 times smaller than the relative tolerance. At
  the default in space, $10^{-4}$, it is about $5\times10^{-6}$ of each
  component's peak, at 25 steps per simulated hour once the start is past;
- after the bulk oxygen falls from 0.21 to 0.05 mol m⁻³, 90 s steps follow
  the change to $3\times10^{-6}$ of air saturation. Six minutes later the
  oxygen profile lies within $7\times10^{-7}$ of the quasi-steady profile on
  the same grid.

Implemented as `marse.core.implicit`, and run by
`marse.core.reactive_transport`.

---

### 9.9 The mouth solved with the box

The pool of §4.8 exchanges with the film on the time scale of the film's
renewal, which is seconds while an intake mixes the film with the mouth. The
S1 prototype first solved the pool apart from the box and corrected it after
each span by exactly what had crossed. That conserved to rounding. But the pH
it gave changed by 0.034 between spans of 60 s and 5 s, and by the same at a
ten times tighter tolerance, so the error was the coupling's. The pool is
therefore solved with the box, in the same implicit steps.

- **Unknowns.** The pool's unknowns are $u = n/H_{\text{ref}}$ for each
  component, where $H_{\text{ref}}$ is the pool's thickness at the resting
  volume. They look like concentrations to the error control and to Newton's
  method, and they are amounts to the update, which books every transfer on
  both sides.
- **Spans.** The mouth runs ahead of the box to the end of each span, at most
  a minute, to its next swallow, or to the next start or end of an intake
  (§4.9). Classical Runge–Kutta in steps of a quarter of a second gives its
  volume, and $H(t)$ is the cubic Hermite polynomial through the samples. A
  swallow ends a span. What was secreted over a step is
  $c_{\text{secreted}}\,[H(t+h) - H(t)]$, so it adds up exactly over a span.
- **A span's own clock.** Each span counts its time from zero. On the run's
  clock, a step of a few seconds taken $t$ hours into a run is the difference
  of two large times, which loses about $\varepsilon\, t/h$ of the step, so
  what the steps add and what the span books drift apart in proportion to
  the time: by $4 \times 10^{-13}$ of the sugar after six hours of meals, and
  more every day. With the span's clock, the box and the mouth balance to
  $3 \times 10^{-16}$ over the same six hours.
- **The flow's lag.** The flow over a span depends on the stimulus in the
  pool, which the box changes as the span goes. The run-ahead expects the box
  to go on giving the stimulus back at the rate it did over the span before.
  Without that, after a rinse, the mouth's sugar changed by 1% between spans
  of 60 s and 5 s, and the plaque's pH by $2 \times 10^{-4}$; with it, by
  $8 \times 10^{-4}$ and $9 \times 10^{-5}$.
- **Exchange as processes.** The exchange with each voxel runs as two
  processes per component: one brings the pool's liquid in at $k\,m$, one
  takes the voxel's back at $k\,c$. The positivity limiter of §9.6 treats the
  second like any consumption. What the pool gives is checked after the
  update, and if the pool would go below zero, it is scaled down and the box
  updated again.
- **The linear systems** are the box's, bordered by the pool:

  $$
  \begin{pmatrix} S & -aB \\ -aC & 1 - aD \end{pmatrix}
  \begin{pmatrix} x_b \\ x_u \end{pmatrix}
  = \begin{pmatrix} r_b \\ r_u \end{pmatrix}.
  $$

  Here $S$ is the box's matrix, $B$ and $C$ couple the pool to the film's
  voxels, and $D$ is the pool's own term, all diagonal in the components. The
  Schur complement on the pool needs $W = S^{-1}B$, one box solve per
  exchanged component, once for each matrix. Every solve after that costs one
  box solve and one $J \times J$ product.

The prototype of this coupling gave the same pH, within $9 \times 10^{-5}$,
with spans of 60 s and 5 s. It conserved the box and the mouth together to
$7 \times 10^{-16}$ over an hour. After a rinse of 10% sucrose held for a
minute, MARSE gives the same pH within $9 \times 10^{-5}$ over an hour with
spans of 60 s and 5 s ([validation](validation.md#the-diet)). Implemented as
`marse.core.reservoir`.

### 9.10 Spreading on voxels

Biomass moves over hours, and solutes relax in milliseconds to seconds. So
the engine alternates: reactions and diffusion run over a spreading interval
$\Delta$ by the implicit method of §9.8, with the biomass growing in place,
and then the excess is moved on. Spreading is a projection back onto
$\phi \le 1$, not a rate, so the splitting is first order in $\Delta$.

**The pressure.** Let $\Omega$ be the voxels with $\phi \ge 1$, and
$e = \phi - 1$ their excess. On $\Omega$, solve the discrete Poisson problem

$$
\sum_{b} (p_a - p_b) = e_a ,
$$

over the face neighbours $b$ of each voxel $a$, with $p = 0$ in every voxel
that has room, and nothing crossing the substratum or the top face. The volume
crossing the face from $a$ to $b$ is $q_{ab} = p_a - p_b$, in voxel volumes, so
each full voxel sends on exactly its excess and what it receives. In a column
full from the substratum, $q$ through each face is the excess summed below
it: the discrete Wanner–Gujer velocity.

**The sweep.** Voxels are visited in order of decreasing pressure. The flow is
the gradient of a potential, so it has no cycles, and in that order each voxel
has received everything flowing into it before it sends. Each moves only
components with a density; the rest stay.

**Order is kept.** A voxel holds a stack of parcels, each a volume and the
amounts in it: its own content, with what flows in from below put beneath it
and what flows in from above put on top. What leaves through the top face is
taken from the top of the stack, and what leaves through the bottom face from
its bottom. A parcel cut in two keeps its amounts in proportion to volume. At
the end, each voxel holds the sum of its parcels.

The simpler choice, mixing each voxel before it sends (donor cell), was built
first and measured. In a film growing without limit, a voxel deep in the film
passes on several times its own volume in each interval, and mixing at every
pass carried a labelled lower layer to the top of the film: after three
doublings, 11% of the material at the top came from the lowest tenth, where
none should have. Its boundary sat 5 to 8 µm low at voxels from 2 µm to
0.5 µm, without converging. With order kept, the boundary sits within a voxel,
and the centroid of the lower layer converges at second order
([validation](validation.md#spreading-in-a-column)). What mixing remains comes
from cutting the voxel that straddles a boundary at each spread, so it grows
with the number of spreads. That is why the interval is not made needlessly
short.

**Rounds.** A voxel with room can receive more than its room. It then joins
$\Omega$, and the pressure and sweep repeat until every voxel fits. Each round
adds voxels to $\Omega$, so the loop ends. A box with no voxel left with room is
refused as full.

**What holds at any interval.**

- *Conservation.* Every transfer leaves one voxel and enters another, so each
  component's total is unchanged to rounding.
- *Positivity.* A voxel never sends more than it holds. A cut takes a
  fraction $f \le 1$ of a parcel, and $x - f x \ge 0$ in floating point too.
  Nothing is clipped.
- *Room.* Every voxel of $\Omega$ ends with exactly what it can hold.

The engine checks these after every spread, whatever mechanism spread:
nothing negative, each component's total unchanged to $10^{-12}$, components
without a density untouched, $\phi \le 1 + 10^{-12}$ everywhere, and the top
layer clear. Then the ledger checks carbon, nitrogen and electrons as after any
step. After a spread, the implicit solver estimates its first step afresh
(§9.8): the moved biomass starts a short transient of its own, and keeping the
previous step had a third of the next steps rejected.

**Measured** in a column ([validation](validation.md#spreading-in-a-column)):
the splitting error in the film's biovolume is first order, with observed
orders 1.05 to 1.22 as $\Delta$ falls from 0.1 h to 0.0125 h. At $\Delta$ = 0.25 h,
the default, a film fed from the liquid thickens within 0.14% of its rate as
$\Delta \to 0$. Implemented as `marse.spatial.spreading`, and run by
`marse.core.reactive_transport`.

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
9. **Carbon, nitrogen, electrons, P, K, Cl and Na only.** Reaction networks
   (§3.6) balance these. Sulphur, calcium, magnesium and other metals are not
   yet balanced, and formulas containing them are refused rather than
   half-checked.
10. **Binding from a stated efficiency.** Cells bind to a surface with an
    attachment efficiency stated per species and material, from a suspension
    that binding does not deplete, and never to cells already bound (§6.4).
11. **pH from local electroneutrality.** The hydrogen ion concentration keeps
    every voxel neutral while ions diffuse independently, with conditional
    constants and no activities (§3.8).
12. **One well-mixed mouth, and a steady diet.** The saliva is one pool,
    secreted and swallowed as Dawes (1983) described, over one site of
    plaque. Drinks and foods enter it steadily over their durations, one at a
    time (§4.8, §4.9).
13. **Spreading by displacement, in a column.** Biomass is incompressible at
    a stated packing density, and moves with one velocity, pushed by the
    excess its growth makes (§6.1, §9.10). Packed plaque settles where it
    shrinks; under the continuum mechanism, growth runs in place for an
    interval before it moves, and a voxel that empties stays empty, since
    collapse and compaction need the matrix and dead biomass of Stage 3.
    Solutes diffuse through biomass as through the liquid. Only packed plaque
    is detached and worn, so a film spreading by the continuum needs a box
    tall enough to hold it.

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
| V4 | Attachment, biofilm initiation | §6.1, §6.3, §6.4 |
| V5 | Two-species competition, $R^{*}$ rule | §7.1, §7.2 |
| V6 | Cross-feeding, coexistence | §3.3, §7.3 |
| V7 | Perturbation and recovery | §2, §8 |
| V8 | Reproducibility from manifest | §9.3, §9.4 |
| S1 | pH, the mouth, the diet and the Stephan curve (criteria G1 to G7) | §3.8, §4.8, §4.9, §9.9 |
| 2d.1 | Oxygen roles (criteria D1 to D3) | §3.9 |
| 2d.2 | Shared space and spreading in a column (criteria D4 to D13) | §6.1, §9.10 |

Full definitions in [`docs/validation.md`](validation.md).

---

## 12. References

Numerical values taken from these sources, with confidence notes, are
tabulated in [`docs/parameters.md`](parameters.md).

1. Alexander, R. (1977) Diagonally implicit Runge–Kutta methods for stiff O.D.E.'s. *SIAM Journal on Numerical Analysis* **14**:1006–1021. [doi:10.1137/0714068](https://doi.org/10.1137/0714068)
2. Alpkvist, E. & Klapper, I. (2007) A multidimensional multispecies continuum model for heterogeneous biofilm development. *Bulletin of Mathematical Biology* **69**:765–789. Bibliographic record not yet checked against the publisher.
3. Andrews, J.F. (1968) A mathematical model for the continuous culture of microorganisms utilizing inhibitory substrates. *Biotechnology and Bioengineering* **10**:707–723. [doi:10.1002/bit.260100602](https://doi.org/10.1002/bit.260100602)
4. Baka, M., Van Derlinden, E., Boons, K., Mertens, L. & Van Impe, J.F. (2013) Impact of pH on the cardinal temperatures of *E. coli* K12: evaluation of the gamma hypothesis. *Food Control* **29**:328–335. [doi:10.1016/j.foodcont.2012.04.022](https://doi.org/10.1016/j.foodcont.2012.04.022)
5. Balaban, N.Q., Merrin, J., Chait, R., Kowalik, L. & Leibler, S. (2004) Bacterial persistence as a phenotypic switch. *Science* **305**:1622–1625. [doi:10.1126/science.1099390](https://doi.org/10.1126/science.1099390)
6. Baranyi, J. & Roberts, T.A. (1994) A dynamic approach to predicting bacterial growth in food. *International Journal of Food Microbiology* **23**:277–294. [doi:10.1016/0168-1605(94)90157-0](https://doi.org/10.1016/0168-1605(94)90157-0)
7. Bardow, A., Moe, D., Nyvad, B. & Nauntofte, B. (2000) The buffer capacity and buffer systems of human whole saliva measured without loss of CO2. *Archives of Oral Biology* **45**:1–12. [doi:10.1016/S0003-9969(99)00119-3](https://doi.org/10.1016/S0003-9969(99)00119-3)
8. Benson, B.B. & Krause, D. (1984) The concentration and isotopic fractionation of oxygen dissolved in freshwater and seawater in equilibrium with the atmosphere. *Limnology and Oceanography* **29**:620–632. [doi:10.4319/lo.1984.29.3.0620](https://doi.org/10.4319/lo.1984.29.3.0620)
9. Bollen, C.M.L., Lambrechts, P. & Quirynen, M. (1997) Comparison of surface roughness of oral hard materials to the threshold surface roughness for bacterial plaque retention: a review of the literature. *Dental Materials* **13**:258–269. [doi:10.1016/S0109-5641(97)80038-3](https://doi.org/10.1016/S0109-5641(97)80038-3)
10. Bolley, C. & Crouzeix, M. (1978) Conservation de la positivité lors de la discrétisation des problèmes d'évolution paraboliques. *RAIRO Analyse numérique* **12**:237–245. [doi:10.1051/m2an/1978120302371](https://doi.org/10.1051/m2an/1978120302371)
11. Briggs, W.L., Henson, V.E. & McCormick, S.F. (2000) *A Multigrid Tutorial*, 2nd edition. SIAM, Philadelphia. [doi:10.1137/1.9780898719505](https://doi.org/10.1137/1.9780898719505)
12. Bruggeman, J., Burchard, H., Kooi, B.W. & Sommeijer, B. (2007) A second-order, unconditionally positive, mass-conserving integration scheme for biochemical systems. *Applied Numerical Mathematics* **57**:36–58. [sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0168927405002242)
13. Busscher, H.J. & van der Mei, H.C. (2006) Microbial adhesion in flow displacement systems. *Clinical Microbiology Reviews* **19**:127–141. [doi:10.1128/CMR.19.1.127-141.2006](https://doi.org/10.1128/CMR.19.1.127-141.2006)
14. Collins, L.M.C. & Dawes, C. (1987) The surface area of the adult human mouth and thickness of the salivary film covering the teeth and oral mucosa. *Journal of Dental Research* **66**:1300–1302. [doi:10.1177/00220345870660080201](https://doi.org/10.1177/00220345870660080201)
15. Crank, J. (1975) *The Mathematics of Diffusion*, 2nd edn. Clarendon Press, Oxford.
16. Dawes, C. & Kubieniec, K. (2004) The effects of prolonged gum chewing on salivary flow rate and composition. *Archives of Oral Biology* **49**:665–669. [doi:10.1016/j.archoralbio.2004.02.007](https://doi.org/10.1016/j.archoralbio.2004.02.007)
17. Dawes, C. & Macpherson, L.M.D. (1992) Effects of nine different chewing-gums and lozenges on salivary flow rate and pH. *Caries Research* **26**:176–182. [doi:10.1159/000261439](https://doi.org/10.1159/000261439)
18. Dibdin, G.H. (1990) Plaque fluid and diffusion: study of the cariogenic challenge by computer modeling. *Journal of Dental Research* **69**:1324–1331. [doi:10.1177/00220345900690062001](https://doi.org/10.1177/00220345900690062001)
19. Dibdin, G.H., Dawes, C. & Macpherson, L.M.D. (1995) Computer modeling of the effects of chewing sugar-free and sucrose-containing gums on the pH changes in dental plaque associated with a cariogenic challenge at different intra-oral sites. *Journal of Dental Research* **74**:1482–1488. [doi:10.1177/00220345950740080801](https://doi.org/10.1177/00220345950740080801)
20. Feder, J. (1980) Random sequential adsorption. *Journal of Theoretical Biology* **87**:237–254. [doi:10.1016/0022-5193(80)90358-6](https://doi.org/10.1016/0022-5193(80)90358-6)
21. Gottlieb, S., Shu, C.-W. & Tadmor, E. (2001) Strong stability-preserving high-order time discretization methods. *SIAM Review* **43**:89–112. [doi:10.1137/S003614450036757X](https://doi.org/10.1137/S003614450036757X)
22. Hairer, E., Nørsett, S.P. & Wanner, G. (1993) *Solving Ordinary Differential Equations I: Nonstiff Problems*, 2nd edition. Springer, Berlin. [doi:10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1)
23. Hairer, E. & Wanner, G. (1996) *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems*, 2nd edition. Springer, Berlin. [doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7)
24. Han, P. & Bartels, D.M. (1996) Temperature dependence of oxygen diffusion in H₂O and D₂O. *Journal of Physical Chemistry* **100**:5597–5602. [doi:10.1021/jp952903y](https://doi.org/10.1021/jp952903y)
25. Heijnen, J.J. & van Dijken, J.P. (1992) In search of a thermodynamic description of biomass yields for the chemotrophic growth of microorganisms. *Biotechnology and Bioengineering* **39**:833–858. [doi:10.1002/bit.260390806](https://doi.org/10.1002/bit.260390806)
26. Henze, M., Gujer, W., Mino, T. & van Loosdrecht, M.C.M. (2000) *Activated Sludge Models ASM1, ASM2, ASM2d and ASM3.* IWA Scientific and Technical Report No. 9. IWA Publishing, London.
27. Hsu, S.-B., Hubbell, S.P. & Waltman, P. (1977) A mathematical theory for single-nutrient competition in continuous cultures of micro-organisms. *SIAM Journal on Applied Mathematics* **32**:366–383. [doi:10.1137/0132030](https://doi.org/10.1137/0132030)
28. Huber, M.L., Perkins, R.A., Laesecke, A. *et al.* (2009) New international formulation for the viscosity of H₂O. *Journal of Physical and Chemical Reference Data* **38**:101–125. [doi:10.1063/1.3088050](https://doi.org/10.1063/1.3088050)
29. Kashket, S., Zhang, J. & Van Houte, J. (1996) Accumulation of fermentable sugars and metabolic acids in food particles that become entrapped on the dentition. *Journal of Dental Research* **75**:1885–1891. [doi:10.1177/00220345960750111101](https://doi.org/10.1177/00220345960750111101)
30. Kovárová-Kovar, K. & Egli, T. (1998) Growth kinetics of suspended microbial cells: from single-substrate-controlled growth to mixed-substrate kinetics. *Microbiology and Molecular Biology Reviews* **62**:646–666. [doi:10.1128/mmbr.62.3.646-666.1998](https://doi.org/10.1128/mmbr.62.3.646-666.1998)
31. Kreft, J.-U., Picioreanu, C., Wimpenny, J.W.T. & van Loosdrecht, M.C.M. (2001) Individual-based modelling of biofilms. *Microbiology* **147**:2897–2912. [doi:10.1099/00221287-147-11-2897](https://doi.org/10.1099/00221287-147-11-2897)
32. Lagerlöf, F. & Dawes, C. (1984) The volume of saliva in the mouth before and after swallowing. *Journal of Dental Research* **63**:618–621. [doi:10.1177/00220345840630050201](https://doi.org/10.1177/00220345840630050201)
33. Lardon, L.A., Merkey, B.V., Martins, S. *et al.* (2011) iDynoMiCS: next-generation individual-based modelling of biofilms. *Environmental Microbiology* **13**:2416–2434. [doi:10.1111/j.1462-2920.2011.02414.x](https://doi.org/10.1111/j.1462-2920.2011.02414.x)
34. Lévêque, A. (1928) Les lois de la transmission de chaleur par convection. *Annales des Mines*, 12th series, **13**:201–299, 305–362, 381–415.
35. Levine, M.J., Reddy, M.S., Tabak, L.A., Loomis, R.E., Bergey, E.J., Jones, P.C., Cohen, R.E., Stinson, M.W. & Al-Hashimi, I. (1987) Structural aspects of salivary glycoproteins. *Journal of Dental Research* **66**:436–441. [doi:10.1177/00220345870660020901](https://doi.org/10.1177/00220345870660020901)
36. Lingström, P. & Birkhed, D. (1993) Plaque pH and oral retention after consumption of starchy snack products at normal and low salivary secretion rate. *Acta Odontologica Scandinavica* **51**:379–388. [doi:10.3109/00016359309040589](https://doi.org/10.3109/00016359309040589)
37. Luedeking, R. & Piret, E.L. (1959) A kinetic study of the lactic acid fermentation. Batch process at controlled pH. *Journal of Biochemical and Microbiological Technology and Engineering* **1**:393–412. [doi:10.1002/jbmte.390010406](https://doi.org/10.1002/jbmte.390010406)
38. Macpherson, L.M.D. & Dawes, C. (1991) Effects of salivary film velocity on pH changes in an artificial plaque containing *Streptococcus oralis*, after exposure to sucrose. *Journal of Dental Research* **70**:1230–1234. [doi:10.1177/00220345910700090101](https://doi.org/10.1177/00220345910700090101)
39. Mei, L., Ren, Y., Busscher, H.J., Chen, Y. & van der Mei, H.C. (2009) Poisson analysis of streptococcal bond-strengthening on saliva-coated enamel. *Journal of Dental Research* **88**:841–845. [doi:10.1177/0022034509342523](https://doi.org/10.1177/0022034509342523)
40. Meinders, J.M., van der Mei, H.C. & Busscher, H.J. (1995) Deposition efficiency and reversibility of bacterial adhesion under flow. *Journal of Colloid and Interface Science* **176**:329–341. [doi:10.1006/jcis.1995.9960](https://doi.org/10.1006/jcis.1995.9960)
41. Monod, J. (1949) The growth of bacterial cultures. *Annual Review of Microbiology* **3**:371–394. [doi:10.1146/annurev.mi.03.100149.002103](https://doi.org/10.1146/annurev.mi.03.100149.002103)
42. Payment, S.A., Liu, B., Offner, G.D., Oppenheim, F.G. & Troxler, R.F. (2000) Immunoquantification of human salivary mucins MG1 and MG2 in stimulated whole saliva: factors influencing mucin levels. *Journal of Dental Research* **79**:1765–1772. [doi:10.1177/00220345000790100601](https://doi.org/10.1177/00220345000790100601)
43. Picioreanu, C., van Loosdrecht, M.C.M. & Heijnen, J.J. (1998) Mathematical modeling of biofilm structure with a hybrid differential-discrete cellular automaton approach. *Biotechnology and Bioengineering* **58**:101–116. [doi:10.1002/(SICI)1097-0290(19980405)58:1<101::AID-BIT11>3.0.CO;2-M](https://doi.org/10.1002/(SICI)1097-0290(19980405)58:1%3C101::AID-BIT11%3E3.0.CO;2-M)
44. Pirt, S.J. (1965) The maintenance energy of bacteria in growing cultures. *Proceedings of the Royal Society B* **163**:224–231. [doi:10.1098/rspb.1965.0069](https://doi.org/10.1098/rspb.1965.0069)
45. Pirt, S.J. (1967) A kinetic study of the mode of growth of surface colonies of bacteria and fungi. *Journal of General Microbiology* **47**:181–197. [doi:10.1099/00221287-47-2-181](https://doi.org/10.1099/00221287-47-2-181)
46. Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (2007) *Numerical Recipes: The Art of Scientific Computing*, 3rd edition. Cambridge University Press, Cambridge. ISBN 978-0-521-88068-8.
47. Prothero, A. & Robinson, A. (1974) On the stability and accuracy of one-step methods for solving stiff systems of ordinary differential equations. *Mathematics of Computation* **28**:145–162. [doi:10.1090/S0025-5718-1974-0331793-2](https://doi.org/10.1090/S0025-5718-1974-0331793-2)
48. Ratkowsky, D.A., Lowry, R.K., McMeekin, T.A., Stokes, A.N. & Chandler, R.E. (1983) Model for bacterial culture growth rate throughout the entire biokinetic temperature range. *Journal of Bacteriology* **154**:1222–1226. [doi:10.1128/jb.154.3.1222-1226.1983](https://doi.org/10.1128/jb.154.3.1222-1226.1983)
49. Ratkowsky, D.A., Olley, J., McMeekin, T.A. & Ball, A. (1982) Relationship between temperature and growth rate of bacterial cultures. *Journal of Bacteriology* **149**:1–5. [doi:10.1128/jb.149.1.1-5.1982](https://doi.org/10.1128/jb.149.1.1-5.1982)
50. Rittmann, B.E. & McCarty, P.L. (2001) *Environmental Biotechnology: Principles and Applications.* McGraw-Hill, New York.
51. Roels, J.A. (1983) *Energetics and Kinetics in Biotechnology.* Elsevier Biomedical Press, Amsterdam.
52. Rosso, L., Lobry, J.R. & Flandrois, J.P. (1993) An unexpected correlation between cardinal temperatures of microbial growth highlighted by a new model. *Journal of Theoretical Biology* **162**:447–463. [doi:10.1006/jtbi.1993.1099](https://doi.org/10.1006/jtbi.1993.1099)
53. Rosso, L., Lobry, J.R., Bajard, S. & Flandrois, J.P. (1995) Convenient model to describe the combined effects of temperature and pH on microbial growth. *Applied and Environmental Microbiology* **61**:610–616. [doi:10.1128/aem.61.2.610-616.1995](https://doi.org/10.1128/aem.61.2.610-616.1995)
54. Saad, Y. & Schultz, M.H. (1986) GMRES: a generalized minimal residual algorithm for solving nonsymmetric linear systems. *SIAM Journal on Scientific and Statistical Computing* **7**:856–869. [doi:10.1137/0907058](https://doi.org/10.1137/0907058)
55. Shellis, R.P. & Dibdin, G.H. (1988) Analysis of the buffering systems in dental plaque. *Journal of Dental Research* **67**:438–446. [doi:10.1177/00220345880670020101](https://doi.org/10.1177/00220345880670020101)
56. Shu, C.-W. & Osher, S. (1988) Efficient implementation of essentially non-oscillatory shock-capturing schemes. *Journal of Computational Physics* **77**:439–471. [doi:10.1016/0021-9991(88)90177-5](https://doi.org/10.1016/0021-9991(88)90177-5)
57. Slot, D.E., Wiggelinkhuizen, L., Rosema, N.A.M. & Van der Weijden, G.A. (2012) The efficacy of manual toothbrushes following a brushing exercise: a systematic review. *International Journal of Dental Hygiene* **10**:187–197. [doi:10.1111/j.1601-5037.2012.00557.x](https://doi.org/10.1111/j.1601-5037.2012.00557.x)
58. Stephan, R.M. (1944) Intra-oral hydrogen-ion concentrations associated with dental caries activity. *Journal of Dental Research* **23**:257–266. [doi:10.1177/00220345440230040401](https://doi.org/10.1177/00220345440230040401)
59. Stewart, P.S. (1998) A review of experimental measurements of effective diffusive permeabilities and effective diffusion coefficients in biofilms. *Biotechnology and Bioengineering* **59**:261–272. [doi:10.1002/(SICI)1097-0290(19980805)59:3<261::AID-BIT1>3.0.CO;2-9](https://doi.org/10.1002/(SICI)1097-0290(19980805)59:3%3C261::AID-BIT1%3E3.0.CO;2-9)
60. Stewart, P.S. (2003) Diffusion in biofilms. *Journal of Bacteriology* **185**:1485–1491. [doi:10.1128/jb.185.5.1485-1491.2003](https://doi.org/10.1128/jb.185.5.1485-1491.2003)
61. Stumm, W. & Morgan, J.J. (1996) *Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters*, 3rd edn. Wiley, New York. ISBN 978-0-471-51185-4.
62. von Ohle, C., Gieseke, A., Nistico, L., Decker, E.M., de Beer, D. & Stoodley, P. (2010) Real-time microsensor measurement of local metabolic activities in ex vivo dental biofilms exposed to sucrose and treated with chlorhexidine. *Applied and Environmental Microbiology* **76**:2326–2334. [doi:10.1128/AEM.02090-09](https://doi.org/10.1128/AEM.02090-09)
63. Walters, M.C., Roe, F., Bugnicourt, A., Franklin, M.J. & Stewart, P.S. (2003) Contributions of antibiotic penetration, oxygen limitation, and low metabolic activity to tolerance of *Pseudomonas aeruginosa* biofilms. *Antimicrobial Agents and Chemotherapy* **47**:317–323. [doi:10.1128/aac.47.1.317-323.2003](https://doi.org/10.1128/aac.47.1.317-323.2003)
64. Wanner, O. & Gujer, W. (1986) A multispecies biofilm model. *Biotechnology and Bioengineering* **28**:314–328. [doi:10.1002/bit.260280304](https://doi.org/10.1002/bit.260280304)
65. Werner, E., Roe, F., Bugnicourt, A. *et al.* (2004) Stratified growth in *Pseudomonas aeruginosa* biofilms. *Applied and Environmental Microbiology* **70**:6188–6196. [doi:10.1128/aem.70.10.6188-6196.2004](https://doi.org/10.1128/aem.70.10.6188-6196.2004)
66. Zwietering, M.H., Jongenburger, I., Rombouts, F.M. & van 't Riet, K. (1990) Modeling of the bacterial growth curve. *Applied and Environmental Microbiology* **56**:1875–1881. [doi:10.1128/aem.56.6.1875-1881.1990](https://doi.org/10.1128/aem.56.6.1875-1881.1990)
67. Zwietering, M.H., Wijtzes, T., de Wit, J.C. & van 't Riet, K. (1992) A decision support system for prediction of the microbial spoilage in foods. *Journal of Food Protection* **55**:973–979. [doi:10.4315/0362-028X-55.12.973](https://doi.org/10.4315/0362-028X-55.12.973)
68. Dawes, C. (1983) A mathematical model of salivary clearance of sugar from the oral cavity. *Caries Research* **17**:321–334. [doi:10.1159/000260684](https://doi.org/10.1159/000260684)
69. Dawes, C. (1989) An analysis of factors influencing diffusion from dental plaque into a moving film of saliva and the implications for caries. *Journal of Dental Research* **68**:1483–1488. [doi:10.1177/00220345890680110301](https://doi.org/10.1177/00220345890680110301)
70. Dawes, C., Watanabe, S., Biglow-Lecomte, P. & Dibdin, G.H. (1989) Estimation of the velocity of the salivary film at some different locations in the mouth. *Journal of Dental Research* **68**:1479–1482. [doi:10.1177/00220345890680110201](https://doi.org/10.1177/00220345890680110201)
71. Dawes, C. & Dibdin, G.H. (1986) A theoretical analysis of the effects of plaque thickness and initial salivary sucrose concentration on diffusion of sucrose into dental plaque and its conversion to acid during salivary clearance. *Journal of Dental Research* **65**:89–94. [doi:10.1177/00220345860650021701](https://doi.org/10.1177/00220345860650021701)
