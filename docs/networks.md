# Reaction networks (configuration schema version 2)

A reaction network lists what a simulation tracks (its **components**) and
how they turn into one another (its **processes**). It is MARSE's
configuration schema version 2, the material core in
[the roadmap](roadmap.md#order-of-work-correctness-first). It is written so
that a process which creates or destroys matter cannot be loaded at all.

With a rate on every process, initial amounts and a clock, a network runs in a
closed, well-mixed box, and every run proves its own balance:

```bash
marse check examples/networks/glucose_cross_feeding.json   # the chemistry, balanced
marse run examples/networks/glucose_cross_feeding.json -o runs/network
marse replay runs/network/manifest.json                    # bit for bit
```

Transport, and with it biofilms, comes next. The version 1 ecosystem
configurations still run with `marse ecosystem` until version 2 replaces them.

## What is checked, and how exactly

Every component declares its chemical formula. From the formula MARSE derives
the three quantities every process must conserve:

| Quantity | Unit | Meaning |
|---|---|---|
| carbon | mol C per mol | carbon atoms |
| nitrogen | mol N per mol | nitrogen atoms |
| electrons | mol e⁻ per mol | the degree of reduction: electrons released on complete oxidation to CO₂, H₂O and NH₄⁺ |

For every process, the carbon, nitrogen and electrons its products hold must
equal what its reactants hold. The check is **exact**. Numbers are read as the
decimals written in the file and the arithmetic is done in fractions, so a
process either balances or it does not; there is no tolerance for a small
leak to hide in. The equations are in
[theory.md §3.6](theory.md#36-composition-continuity-and-the-degree-of-reduction).

Water and protons are never listed. Once carbon, nitrogen and electrons
balance, the oxygen, hydrogen and charge balances can always be closed with
H₂O and H⁺, and MARSE does so implicitly. `marse check` shows the amounts.

## A first network

`examples/networks/glucose_cross_feeding.json` describes two bacteria that
share glucose. An aerobic heterotroph respires it. A lactic fermenter ferments
it, and the heterotroph also grows on the fermenter's lactate. Its first
growth process:

```json
{
  "name": "heterotroph_growth_on_glucose",
  "kind": "growth",
  "biomass": "heterotroph",
  "substrate": "glucose",
  "yield_mol_per_mol": 3.6,
  "balanced_by": ["oxygen", "carbon_dioxide", "ammonium"]
}
```

It states only the yield: 3.6 mol of biomass (C-mol, one carbon each) per mol
of glucose. The balances fix everything else. `marse check` prints

```text
heterotroph_growth_on_glucose (growth of heterotroph on glucose, per mol formed)
  0.277778 glucose + 0.616667 oxygen + 0.2 ammonium -> 0.666667 carbon_dioxide +
    heterotroph + 1.06667 H2O + 0.2 H+
```

The fermenter's lactate is derived in the same way. Its process lists
`lactate` in `balanced_by`, so the amount of lactate is whatever carbon and
electrons the glucose supplies beyond the new biomass. A product that no
substrate pays for, the fault in the version 1 engine
([known defect 2](validation.md#the-two-dimensional-ecosystem-engine-is-not-yet-verified)),
cannot be written down.

## The file

A network is one JSON object. Numeric field names end in their unit (see
[Units](#units-are-part-of-the-names)), and every key is checked: an unknown
key is an error, not a silently ignored typo.

### Network fields

| Field | Type | Required | Meaning |
|---|---|---|---|
| `schema_version` | integer | yes | Always `2`. A file without it is taken to be version 1 and pointed at `marse ecosystem`. |
| `description` | text | no | What the network represents, and where its numbers come from. |
| `components` | list of components | yes | At least one. |
| `processes` | list of processes | yes | May be empty. |
| `experiment_id` | text | to run | Names the run in its manifest. |
| `initial_mol_per_m3` | object of numbers | no, default `0` each | The starting concentration of each component, in mol per m³ (mmol per litre). Components left out start at zero. |
| `duration_h` | number | to run | How long to run. |
| `timestep_h` | number | to run | The longest substep, and the unit of recording. Accuracy comes from the tolerances, not from this. |
| `record_interval_h` | number | no, default `timestep_h` | How often the trajectory records a row, rounded to a whole number of timesteps. The final state is always recorded. |
| `relative_tolerance` | number | no, default `1e-6`, or `1e-4` in space | The accepted local error, as a fraction of each component's largest value. In space that is its largest value anywhere in the box. |
| `absolute_tolerance_mol_per_m3` | number | no, default `1e-9` | The accepted local error for traces (picomolar). |
| `seed` | integer | no, default `0` | Places `random_colonies`, and is recorded in the manifest. Nothing else draws random numbers. |
| `domain` | domain | no | Runs the network in space instead of a closed box; see [Running a network in space](#running-a-network-in-space). |

### Component fields

| Field | Type | Required | Meaning |
|---|---|---|---|
| `name` | name | yes | Letters, digits and underscores, starting with a letter, at most 64 characters. Names become keys and file names in outputs. |
| `phase` | `dissolved` or `particulate` | yes | Dissolved components are carried by the liquid. Particulate ones, such as biomass, move only with the biofilm. |
| `formula` | text | yes | The chemical formula of one unit of the component, such as `C6H12O6` or `CH1.8O0.5N0.2`. |
| `charge` | number | no, default `0` | The charge of one formula unit, in elementary charges: `1` for ammonium, `-1` for lactate. |

A component is counted in mol of its formula unit. Biomass written per carbon
atom, as `CH1.8O0.5N0.2`, is therefore counted in C-mol. The formula rules:

- **Elements** may be C, H, O and N, each followed by an optional count. A
  count may be a decimal, for a lumped composition such as biomass. A symbol
  may repeat (`CH3COOH`) and its counts are summed. Sulphur and phosphorus
  arrive in Stage 3; parentheses are not supported.
- **The charge goes in its own field.** A formula whose counts are whole
  numbers must have an even number of electrons. An odd number almost always
  means an ion written without its charge: `NH4` without `"charge": 1` would
  give ammonium an electron it does not have, and every balance through it
  would then be wrong. Radicals, which really have an odd number, are not
  supported yet.
- **Water, protons and hydroxide are refused.** They hold no carbon, nitrogen
  or electrons, so no balance could constrain them; they are closed implicitly
  instead.

### Growth process fields

| Field | Type | Required | Meaning |
|---|---|---|---|
| `name` | name | yes | Unique among processes. |
| `kind` | `growth` | yes | |
| `biomass` | component name | yes | What grows; must be particulate. |
| `substrate` | component name | yes | What it grows on. |
| `yield_mol_per_mol` | number | yes | mol of biomass formed per mol of substrate consumed; positive. |
| `products_mol_per_mol` | object of numbers | no | mol of each product released per mol of substrate consumed, for products whose yield is known; positive. |
| `balanced_by` | list of component names | no | Components whose coefficients the balances determine. |
| `rate` | rate | to run | How fast the process runs; see [Rates](#rates). |

A growth process is written per mol of biomass formed: biomass +1, substrate
−1/Y, each product +Y<sub>P</sub>/Y, and the `balanced_by` components solved.
So the rate of the process is the growth rate of the biomass.

### Reaction process fields

| Field | Type | Required | Meaning |
|---|---|---|---|
| `name` | name | yes | Unique among processes. |
| `kind` | `reaction` | yes | |
| `stoichiometry_mol_per_mol` | object of numbers | yes | The coefficient of each component per unit of reaction: negative if consumed, positive if produced, never zero. |
| `balanced_by` | list of component names | no | Components whose coefficients the balances determine. |
| `rate` | rate | to run | How fast the process runs; see [Rates](#rates). |

Endogenous respiration of the heterotroph, for example, gives only the biomass
consumed and lets the balances find the oxygen, carbon dioxide and ammonium:

```json
{
  "name": "heterotroph_endogenous_respiration",
  "kind": "reaction",
  "stoichiometry_mol_per_mol": {"heterotroph": -1},
  "balanced_by": ["oxygen", "carbon_dioxide", "ammonium"]
}
```

## Rates

A process runs at

$$
r = k \; c_{\text{proportional\_to}} \prod_i f_i ,
$$

in mol of the process's reference per m³ per hour: for growth, mol of biomass
formed. Each factor f is dimensionless. The equations are in
[theory.md §3.7](theory.md#37-rates).

### Rate fields

| Field | Type | Required | Meaning |
|---|---|---|---|
| `maximum_per_h` | number | yes | k, the rate per unit of `proportional_to` when every factor is 1; not negative. |
| `proportional_to` | component name | growth: no, default the biomass; reaction: yes | The component the rate is proportional to. |
| `factors` | list of factors | no | The dimensionless terms that speed or slow the process. |
| `assumed_in_excess` | list of component names | no | Components the process consumes although its rate does not depend on them. |

### Factor fields

| Field | Type | Required | Meaning |
|---|---|---|---|
| `component` | component name | yes | The component the factor responds to. |
| `form` | `monod`, `inhibition` or `haldane` | yes | `monod`: S/(K+S). `inhibition`: K_I/(K_I+S). `haldane`: S/(K+S+S²/K_I). |
| `half_saturation_mol_per_m3` | number | monod, haldane | K, positive. |
| `inhibition_mol_per_m3` | number | inhibition, haldane | K_I, positive. |

Two rules are checked when the file is read:

- **Each component limits a process at most once.** A second factor for the
  same component is refused. This is the error behind known defect 4, where
  the version 1 engine applied a substrate's Monod term twice. A substrate
  that also inhibits takes one `haldane` factor.
- **A process may only consume what its rate depends on.** Every component a
  process consumes must be its `proportional_to`, or carry a `monod` or
  `haldane` factor, so that the process slows to a stop as the component runs
  out. If a component really never runs short, list it in
  `assumed_in_excess`; activated-sludge models make that assumption for
  ammonium, for example. The assumption is then recorded, and a run in which
  it fails stops with an error rather than continuing on a false premise.

## Running a network

`marse run` integrates the network in a closed, well-mixed box, from
`initial_mol_per_m3` for `duration_h`. It writes `trajectory.csv`, with units
in the headers, and `manifest.json`. On the example:

```text
kind        well_mixed
steps       96 over 48 h: 7952 adaptive substeps, 9 retried, 53 limited to keep concentrations positive
final       mol per m3
  glucose                  1.37689e-21
  oxygen                   1.16245e-275
  ammonium                 7.54085
  carbon_dioxide           0.864787
  lactate                  35.6132
  heterotroph              0.34708
  fermenter                11.9687
balance     carbon, nitrogen and electrons conserved to 8.3e-15 of their totals
```

In the closed box, the heterotroph uses up the oxygen within eight hours.
The fermenter then turns the remaining glucose into lactate. Cross-feeding on
that lactate needs a continuing oxygen supply, which transport will bring.

Three guarantees hold for every run:

- **Matter is conserved.** After every step a ledger compares the carbon,
  nitrogen and electron totals with the start. A drift beyond 1e-9 of the
  total stops the run with an error. The largest drift seen is written to the
  manifest, and is typically around 1e-14.
- **Concentrations stay non-negative, without clipping.** When a step would
  take more of a component than there is, only the processes consuming it are
  slowed, each as a whole, so their balances are untouched.
- **The run replays bit for bit.** The substeps are chosen by a
  deterministic error control. `marse replay` recomputes the run and compares
  a digest of the final state.

The method is in [theory.md §9.6](theory.md#96-positive-conservative-integration).

## Running a network in space

With a `domain`, the same network runs in a box of cubic voxels over a flat
surface, in one, two or three dimensions. The scene is the standard one for a
biofilm:

- **Colonies sit on the surface.** The surface might be a tooth, a pipe wall or
  a slide.
- **Liquid lies above them.** Its bulk concentrations are held at the top face
  of the box.
- **The box repeats sideways.** Its lateral faces are periodic, so it is one
  tile of a larger surface.

Dissolved components diffuse between voxels at the stated diffusivities. Every
process runs in every voxel at once, at the rates described in [Rates](#rates). Biomass
grows where it is. Colonies spreading, and sharing space as they do, come in the
next stage.

```bash
marse check examples/networks/surface_biofilm_3d.json   # the space, and what an explicit step would cost
marse run examples/networks/surface_biofilm_3d.json -o runs/surface
marse replay runs/surface/manifest.json
```

`marse run` writes:

- `totals.csv`: each component per m² of surface over time, and what has
  entered through the top;
- `manifest.json`;
- `frames/`: every component's field at each recorded time, as NumPy files;
- `vtk/`: the same frames for 3-D viewers. Open `vtk/run.pvd` in
  [ParaView](https://www.paraview.org/) and step through time. A contour of
  oxygen shows where colonies have made the liquid anoxic.

### The example

`surface_biofilm_3d.json` places a heterotroph colony and a fermenter colony,
each 40 µm in radius, on a 160 µm square of surface. Their centres are 60 µm
apart, so the two overlap near the surface. Six small heterotroph colonies are
scattered from the seed. The liquid above holds a quarter
of air-saturated oxygen. `marse run` prints:

```text
experiment  surface-biofilm-3d
kind        reactive_transport
space       3-D, 32 x 32 x 16 voxels of 5 um (160 x 160 x 80 um)
steps       8 over 2 h: 343 implicit substeps, 1 retried, 78 limited to keep concentrations positive
final                      mol per m2 of surface   entered through the top
  glucose                   0.000364063       +0.01506
  oxygen                    2.23469e-06      +0.004527
  ammonium                  0.000399011      +0.002861
  carbon_dioxide            8.21439e-05      -0.005242
  lactate                   3.41743e-05       -0.02363
  heterotroph                 0.0173451             +0
  fermenter                   0.0190625             +0
balance     carbon, nitrogen and electrons, counting what crossed the top, conserved to 1.8e-16
```

After two hours:

- **The heterotroph colony has an anoxic core.** Oxygen at its base is 0.03%
  of the liquid's, and 375 voxels hold less than 1%.
- **The fermenter makes lactate**: up to 1.5 mol m⁻³ at its base. Its biomass
  has grown 1.8-fold in place, the heterotrophs' 1.45-fold.
- **Lactate feeds the heterotrophs.** It diffuses to them from the fermenter,
  and 39% of their growth now runs on it. Most still escapes to the liquid,
  0.024 mol m⁻² over the two hours. At the end the heterotrophs take up
  lactate at about a twentieth of the rate it leaves.

The run takes about four and a half minutes on one core. Most of that is the
first quarter hour, in 302 of the 343 steps. The colonies start in fresh
liquid, use up the oxygen within seconds, and lactate appears from nothing: a
transient spread over five decades of time ([theory.md
§9.8](theory.md#98-implicit-reactiontransport-integration)). After it, steps
grow to about three minutes. `surface_biofilm_1d.json` runs the same
chemistry in a single column in three seconds.

### Domain fields

| Field | Type | Required | Meaning |
|---|---|---|---|
| `voxels` | list of 1 to 3 whole numbers | yes | Voxels along each axis. The last axis is height above the surface. `[64]` is a column, `[32, 64]` a vertical slice, `[32, 32, 16]` a box. Counts divisible by 2 several times make the solver fastest. |
| `voxel_um` | number | yes | The edge of a cubic voxel, in µm. |
| `bulk_mol_per_m3` | object of numbers | no, default `0` each | What the liquid above holds, held at the top face. Dissolved components only. |
| `diffusivity_m2_per_s` | object of numbers | yes | One value for every dissolved component, in m² per s, as tables give them. Particulate components (biomass) do not diffuse. |
| `colonies` | list of colonies | no | Hemispheres of biomass on the surface, placed where stated. |
| `random_colonies` | list of random colonies | no | Hemispheres of biomass scattered over the surface from the run's `seed`. |

`initial_mol_per_m3` fills every voxel. Each colony then sets its component to
its concentration in the voxels whose centres lie inside it. A colony that
covers no voxel is refused.

### Colony fields

| Field | Type | Required | Meaning |
|---|---|---|---|
| `component` | component name | yes | A particulate component, such as a species' biomass. |
| `center_um` | list of numbers | yes | Where it stands on the surface, one coordinate per lateral axis, in µm: two in 3-D, one in 2-D, none in 1-D, where a colony is a layer. |
| `radius_um` | number | yes | Its radius, which is also its height. |
| `concentration_mol_per_m3` | number | yes | The biomass density inside it, in C-mol per m³ for biomass written per carbon atom. |

### Random colony fields

| Field | Type | Required | Meaning |
|---|---|---|---|
| `component` | component name | yes | A particulate component. |
| `count` | integer | yes | How many to scatter. |
| `radius_um` | number | yes | Each one's radius. It must be large enough to always cover a voxel. |
| `concentration_mol_per_m3` | number | yes | The density inside each. |

The positions are drawn from the random stream `colonies`, derived from `seed`
and recorded in the manifest. So a replay puts every colony back in the same
voxels, and another seed scatters them anew.

Every run in space keeps the three guarantees of a closed box:

- **Matter is conserved, counting what crosses the top.** After every step, a
  ledger compares the carbon, nitrogen and electrons in the box with what has
  entered or left through the top face. That exchange is computed from the face
  fluxes themselves, not inferred from the difference. A drift beyond 1e-9
  stops the run.
- **Concentrations stay non-negative, without clipping.**
- **The run replays bit for bit.**

The method is implicit, so a step can be minutes long even though diffusion
across a voxel takes milliseconds. It follows fast changes, such as a sudden
drop in the oxygen above, and settles onto the quasi-steady profile when
nothing is changing fast ([theory.md §9.8](theory.md#98-implicit-reactiontransport-integration)).

## Balancing: `balanced_by`

Each component in `balanced_by` gets the coefficient that makes the process
balance. There are three balances (carbon, nitrogen and electrons), so at most
three components can be determined, and they must be determined uniquely.
MARSE refuses:

- **a quantity nothing in `balanced_by` carries.** Growth on glucose balanced
  only by oxygen and carbon dioxide leaves the biomass's nitrogen with no
  source. The message names nitrogen and suggests the nitrogen source.
- **components the balances cannot tell apart.** Glucose and lactate both hold
  four electrons per carbon and no nitrogen, so balancing with both is
  ambiguous.
- **a set that cannot balance everything at once**, such as a single biomass
  asked to absorb glucose's carbon and electrons in the wrong ratio.
- **a component that is both given and in `balanced_by`.**

A process without `balanced_by` must balance exactly as written. A copied
coefficient rounded to 5.999 instead of 6 is refused. Let MARSE derive it
instead of typing it.

A component in `balanced_by` may come out as zero. In homolactic fermentation,
glucose balanced by lactate and carbon dioxide gives 2 lactate and no carbon
dioxide, which `marse check` reports.

## Worked example: aerobic growth on glucose

Biomass CH<sub>1.8</sub>O<sub>0.5</sub>N<sub>0.2</sub> holds 4 + 1.8 − 2(0.5) −
3(0.2) = 4.2 electrons per C-mol; glucose holds 24 per mol; O₂ holds −4. With
a yield of 3.6 C-mol per mol of glucose, per C-mol of biomass formed:

| Balance | Equation | Solved |
|---|---|---|
| — | glucose = −1/3.6 | −5/18 mol |
| carbon | −6 × 5/18 + 1 + ν<sub>CO₂</sub> = 0 | ν<sub>CO₂</sub> = +2/3 |
| nitrogen | 0.2 + ν<sub>NH₄⁺</sub> = 0 | ν<sub>NH₄⁺</sub> = −1/5 |
| electrons | −24 × 5/18 + 4.2 − 4 ν<sub>O₂</sub> = 0 | ν<sub>O₂</sub> = −37/60 |

Oxygen atoms then give 16/15 mol of water produced, and charge gives 1/5 mol
of H⁺ released:

5/18 C₆H₁₂O₆ + 37/60 O₂ + 1/5 NH₄⁺ → CH₁.₈O₀.₅N₀.₂ + 2/3 CO₂ + 16/15 H₂O + 1/5 H⁺

Of the 6.67 mol of electrons in the glucose, 4.2 (63%) end up in biomass and
the rest go to oxygen. The test suite checks this row exactly. It also checks
growth rows against the half-reaction method of Rittmann and McCarty (2001,
chapter 2), an independent route to the same stoichiometry, and it checks
textbook reactions: respiration, alcoholic and homolactic fermentation, and
nitrification.

## Units are part of the names

Every numeric field ends in its unit, so a value can never be read in the
wrong one. A key with a different unit, such as `yield_g_per_g`, is pointed at
the field that exists.

| Suffix | Unit |
|---|---|
| `_mol_per_mol` | mol of one component per mol of another |
| `_mol_per_m3` | mol per cubic metre, which is mmol per litre |
| `_m2_per_s` | square metres per second, as diffusivities are tabulated |
| `_per_h` | per hour |
| `_um` | micrometres |
| `_h` | hours |

Six numeric fields have no unit: `schema_version`, a version number;
`charge`, in elementary charges; `seed`, an identifier; `relative_tolerance`,
a fraction; and `voxels` and `count`, which are counts. The test suite checks that every numeric
field follows this rule, and that the tables on this page list exactly the
fields MARSE reads.

## Converting a literature yield

Yields are often published in grams of dry biomass per gram of substrate.
Convert with the molar masses, which `marse check` prints:

$$
Y_{\mathrm{mol/mol}} = Y_{\mathrm{g/g}} \times \frac{M_{\mathrm{substrate}}}{M_{\mathrm{biomass}}},
\qquad \text{e.g. } 0.5 \times \frac{180.156}{24.626} = 3.66 .
$$

Dry weight includes ash, which the biomass formula does not. If the published
yield counts ash, multiply by the ash-free fraction first. Record where every
yield comes from in the network's `description`, with its confidence level
from [parameters.md](parameters.md#how-to-read-this-table). The version 2
schema will carry these grades itself once the examples are rebuilt on it.

## What is not in a network yet

The next increments, in the order of
[the roadmap](roadmap.md#order-of-work-correctness-first):

- biomass that spreads as it grows and shares space between species, with
  oxygen roles for species and heritable lineages;
- the examples and the periodontal study rebuilt on version 2, after which
  version 1 is removed;
- dead biomass, the extracellular matrix, and a diffusivity that depends on
  them.
