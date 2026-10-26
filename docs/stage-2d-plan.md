# Stage 2d plan: biomass that spreads and shares space

**Status: in progress.** 2d.1, oxygen roles, is done: it met D1 to D3
([validation.md](validation.md#oxygen-roles)). 2d.2, shared space and
spreading in a column, is done: it met D4 to D13
([validation.md](validation.md#spreading-in-a-column)), with the changes
recorded under its criteria. Nothing else is built yet. The
criteria below are set before the code is written, as they were for 2c and
S1, so that the stage can fail them. When an increment lands, its results
replace its thresholds in [validation.md](validation.md), and this page
records what changed from the plan.

Stage 2d is the fourth increment of the material core
([roadmap](roadmap.md#order-of-work-correctness-first)). Its deliverables are
biomass that spreads and shares space in 3-D, an oxygen role for each species,
and variants as heritable lineages. The roadmap also asks this plan to choose
between a continuum and individual cells.

## Why 2d comes next

- **Colonies grow in place today.** The 64 × 64 × 32 performance run ended
  with about 55 times its starting biomass in the same voxels
  ([validation.md](validation.md#performance)). Nothing biological can be
  read from a run in space until biomass moves.
- **It holds the requirements of four known defects.** The engine already
  meets defects 2, 4 and 7, and 2e meets defect 1 when it rebuilds the
  examples. 2d meets the other four, listed below.
- **Every later stage needs it.** 2e rebuilds the examples on spreading
  biofilms. S2's five species compete for space. Stage 5's IWA benchmarks
  BM1 and BM3 are films that thicken. Stage 6's periodontal study depends on
  where each species ends up.

| # | Defect in the version 1 engine | Requirement on the version 2 engine | Increment |
|---|---|---|---|
| 3 | The periodontal anaerobes need oxygen to grow | A species declared anaerobic cannot depend on oxygen, and grows without it | 2d.1 |
| 5 | Moving biomass clips negative values, which creates matter | Spreading is conservative and positive by construction, at any interval | 2d.2 |
| 8 | Each species has its own carrying capacity | No voxel holds more biomass than fits, all species together | 2d.2 |
| 6 | Mutation marks grid cells, not lineages | A variant exists only where cells of it are, and moves with them | 2d.4 |

The expected-failure tests in `tests/test_ecosystem_known_defects.py` run the
version 1 engine, so they stay as they are until 2e removes it. 2d writes each
requirement as an ordinary test on the version 2 engine.

## The decision: a continuum first, individual cells later

[modeling-landscape.md §2](modeling-landscape.md#2-the-biomass-spreading-decision)
found that the spreading mechanism changes conclusions about competition and
cooperation. It therefore recommends making spreading a declared, swappable
provider, and reporting a sensitivity check across at least two mechanisms.
This plan follows that. The choice left is **which mechanisms come first**.

**Recommendation: two mechanisms on the voxel grid in 2d, a continuum and a
cellular automaton. Individual cells come later, as a third provider behind
the same contract.**

The reasons:

- **Every guarantee of the version 2 engine is stated per voxel.** That
  covers the ledger, positivity without clipping, bit-for-bit replay and the
  pH. A voxel mechanism moves biomass by transfers between neighbouring
  voxels, the same structure as diffusion, so every guarantee carries over
  unchanged. Individual cells would need a second representation, mapped to
  voxels and back every step, before any of them held again.
- **A continuum has a one-dimensional limit, and individual cells have
  none.** In a column, the continuum mechanism is exactly the displacement
  model of Wanner and Gujer (1986). That model is the reference for
  everything in Stages S and 5:
  - S4 runs years of daily cycles in 1-D columns.
  - All nine published solutions to IWA BM3 are one-dimensional.
  - BM1 is the 1-D V3 case with a film that grows.
- **The scenes already exist on voxels.** E1's bound cells are C-mol in the
  bottom voxels, and S1's plaque is a column of voxels. A continuum spreads
  them without changing either.
- **Cost.** The 2c engine already misses its performance target by a factor
  of 2.3. Mapping agents onto voxels would add cost per step on top of that.
  A spreading step on the grid costs about as much as one pressure solve.

What this gives up, stated so that it is not forgotten:

- **Cell shape, and mechanics with attraction.** The literature ranks
  force-based individual cells as the mechanism that mixes lineages least.
  The continuum here is the voxel analogue of mechanical relaxation, not a
  replacement for it.
- **Mutation as a rare, random event.** In a continuum, a mutant appears as
  a continuous flux (2d.4). The random arrival of the first mutant cell waits
  for individual cells.
- **The sensitivity check has a narrow span.** Continuum against cellular
  automaton spans less of the range than continuum against individual cells
  would. Any spatial claim made before the third provider exists must say
  so.

**This is the decision the project's owner needs to confirm.** The plan below
assumes it. If individual cells are wanted in v1.0, the provider contract of
2d.2 still applies, but an increment 2d.6 has to be added and the stage
grows by about half.

## The model

### Space is shared: volume fractions

Every particulate component that takes up room declares a **packing density**
ρ<sub>j</sub>, in mol per m³. That is its concentration when it alone fills a
voxel. A voxel's biomass then fills the fraction

$$
\phi = \sum_j \frac{c_j}{\rho_j} \le 1 ,
$$

summed over every particulate component with a density. The capacity is shared
by construction: two species together cannot exceed it, which is the
requirement behind defect 8. After every spreading step, the engine checks
φ ≤ 1 in every voxel, alongside the ledger. A violation stops the run, as a
leak does.

- **Initial states.** An initial state that overfills a voxel is refused
  when the file loads.
- **Components with no density take no room and do not move.** They include
  food retained on the teeth, which lies in the film, and the reversibly
  bound cells of E1, which are locked or released within minutes. A test
  fixes each case.
- **What a reference density looks like.** IWA BM1 gives its film a density
  of 10 kg COD per m³. Biomass CH<sub>1.8</sub>O<sub>0.5</sub>N<sub>0.2</sub>
  holds 4.2 e⁻ per C-mol, which is 33.6 g COD per C-mol, so that is about
  300 C-mol per m³. This figure has to be checked against the primary text
  before it is graded in parameters.md.

### Spreading: moving the excess

Growth runs in place during a span of the implicit solver, exactly as it does
today. Afterwards, some voxels hold more than they can. Spreading then moves
the excess until every voxel fits again. For the continuum mechanism:

1. **The full region.** Let Ω be the voxels with φ ≥ 1, and e = φ − 1 the
   excess in each of them.
2. **Pressure.** Solve the discrete Poisson problem −∇<sub>h</sub>²p = e on Ω.
   The boundary conditions:
   - p = 0 in every voxel that has room;
   - no flux through the substratum;
   - periodic lateral faces, as for diffusion.

   The volume moving across the face from voxel a to voxel b is
   q<sub>ab</sub> = p<sub>a</sub> − p<sub>b</sub>, in voxel volumes. Biomass
   flows from high pressure to low, as an incompressible material through a
   porous medium does (Darcy's law). This is the multidimensional continuum
   model of Alpkvist and Klapper (2007), written on voxels.
3. **The sweep.** Visit the voxels of Ω in order of decreasing pressure. Each
   voxel first receives what flows in from higher pressure. It then sends
   q<sub>ab</sub> to each lower neighbour, carrying every component in
   proportion to what the voxel now holds. The flow is the gradient of a
   potential, so it has no cycles, and decreasing pressure is an order in
   which every voxel sends only after it has received everything.
4. **Repeat.** A voxel with room may receive more than its room. If so, it
   joins Ω and steps 1 to 3 repeat. Each round adds voxels to Ω, so the loop
   ends, usually after one to three rounds.

The sweep's guarantees hold at **any** interval, with no stability limit:

- **Conservation.** Every transfer leaves one voxel and enters its
  neighbour, and particulates never cross the top face or the substratum. The
  ledger therefore holds whatever the transfers are.
- **Positivity.** A voxel of Ω holds 1 + e + inflow and sends e + inflow,
  which is never more than it holds.
- **Single velocity.** All the components in a voxel move together. The
  composition a voxel sends is the one it holds, so a species, or a variant,
  appears only where its cells came from (2d.4).

**In one dimension** the Poisson problem integrates directly: the flux
through each face is the excess summed below it. That is the Wanner–Gujer
displacement velocity, u(z) = ∫₀<sup>z</sup> Σ<sub>j</sub> r<sub>j</sub>/ρ<sub>j</sub> dz′,
in discrete form. A test holds the two to rounding (criterion D6).

**A voxel that empties stays empty.** Decay or lysis that lowers φ leaves a
void, and the biomass around it does not move in to fill it. Collapse and
compaction need the matrix and dead biomass of Stage 3. Until then, this is
an assumption listed in theory.md §10.

**The top of the box.** Biomass that would spread into the top layer stops
the run with an error, which names the time and asks for a taller box.
Detachment, which would balance growth, is Stage 3.

### Coupling to the implicit solver

Solutes relax in milliseconds to seconds, and biomass moves over hours. The
engine therefore alternates:

1. integrate reactions and diffusion over a spreading interval with the
   implicit method of [theory.md §9.8](theory.md#98-implicit-reactiontransport-integration),
   unchanged;
2. spread;
3. check the ledger and the capacity.

Spreading is a projection back onto φ ≤ 1, not a rate, so this splitting is
first order in the interval. Biomass grows for one interval where it was
before it moves. Criterion D7 measures that error, and the default interval
is chosen from the measurement. The solver's adaptive substeps, its error
control and its limiter are not touched.

### The provider contract

A spreading provider receives the state after a span and returns face
transfers of particulate components. The engine, not the provider, checks
four things after every step:

- the ledger;
- φ ≤ 1 in every voxel;
- no negative value;
- no transfer through the top face or the substratum.

So a new mechanism cannot weaken a guarantee by mistake. The manifest
records the provider's name and version (`spreading: continuum_pressure_v1`).
A run without a `spreading` block grows in place exactly as now, and its
results do not change by a single bit.

### The second mechanism: a cellular automaton

This follows Picioreanu, van Loosdrecht and Heijnen (1998):

- **What moves.** Each over-full voxel moves its excess along a shortest
  path, through face neighbours, to the nearest voxel with room.
- **How.** Each voxel on the path passes the same volume on to the next,
  carrying its own composition.
- **Order and ties.** Over-full voxels are taken in order of decreasing
  excess, then by index. Paths of equal length are chosen from the seeded
  stream `spreading`, which is recorded in the manifest.

Conservation and positivity follow by the same argument as for the sweep.
This mechanism shoves along one path where the continuum spreads in every
direction, so it should mix lineages more. Criterion D26 measures how much
more.

### Oxygen roles

Every species, meaning every component that is the `biomass` of a growth
process, declares an `oxygen_role`. Oxygen is the dissolved component whose
formula is O<sub>2</sub>, so no new field is needed to name it. Each role is
a rule the network must obey when it loads:

| Role | Its processes |
|---|---|
| `obligate_aerobe` | Every growth process consumes oxygen. |
| `microaerophile` | Every growth process consumes oxygen and has a `haldane` factor on it, so that too much oxygen slows it. |
| `facultative` | At least one growth process consumes oxygen and at least one does not. |
| `aerotolerant` | No process consumes oxygen, and none has a `monod` or `haldane` factor on it. |
| `obligate_anaerobe` | As `aerotolerant`, and every growth process has an `inhibition` factor on oxygen. |

- **Defect 3 cannot be written down.** The version 1 periodontal anaerobes
  had a Monod term on oxygen, and that is refused for any anaerobe. A refusal
  names the process and the factor.
- **"Its processes".** A species' processes are its growth processes and the
  reactions that are `proportional_to` it.
- **Variants inherit the role** of the species they come from.

### Variants and lineages

A **variant** is a particulate component declared `variant_of` another, with
the same formula:

- **Its processes are copies of the parent's.** It gets every growth process
  and every reaction `proportional_to` the parent, with the biomass renamed.
  The copies are named after the parent's, with the variant's name appended.
- **What it may change.** A `processes` map, keyed by the parent's process
  names, may change rate fields and the yield. Since the formula is the same,
  each copy is balanced again when it loads.
- **A variant with no changes is neutral.** It is a label. Neutral variants
  seeded as separate colonies are how 2d.4 and 2d.5 measure lineage mixing.

A **mutation** is a new kind of process, from a parent to a variant:

```json
{"name": "heterotroph_to_resistant", "kind": "mutation",
 "from": "heterotroph", "to": "heterotroph_resistant", "per_division": 1e-6}
```

- **Its rate** is `per_division` times the summed rate of the parent's growth
  processes. It is balanced because both components have the same formula,
  and its Jacobian is the same multiple of theirs.
- **Where it happens.** The rate is proportional to the parent's growth, so a
  mutation happens only where the parent is dividing. It is exactly zero
  where the parent is absent, which is the requirement behind defect 6.
- **It moves with its cells.** Spreading carries a variant with the biomass
  around it.
- **Closed form.** For a neutral variant in a closed box, with growth rate μ
  and mutation fraction p, the parent follows P₀e<sup>(1−p)μt</sup> and the
  variant P₀(e<sup>μt</sup> − e<sup>(1−p)μt</sup>). Criterion D21 tests
  against this.

## The configuration

New fields, each with its unit in its name, following
[networks.md](networks.md#units-are-part-of-the-names):

```json
{
  "components": [
    {"name": "heterotroph", "phase": "particulate", "formula": "CH1.8O0.5N0.2",
     "density_mol_per_m3": 300, "oxygen_role": "facultative"},
    {"name": "heterotroph_resistant", "phase": "particulate", "formula": "CH1.8O0.5N0.2",
     "density_mol_per_m3": 300, "variant_of": "heterotroph",
     "processes": {"heterotroph_growth_on_glucose": {"maximum_per_h": 0.25}}}
  ],
  "domain": {
    "spreading": {"mechanism": "continuum", "interval_h": 0.05}
  }
}
```

| Where | Field | Meaning |
|---|---|---|
| component | `density_mol_per_m3` | The packing density ρ. Required on every particulate component that moves with the biofilm, when the domain spreads. |
| component | `oxygen_role` | One of the five roles. Required on every species. |
| component | `variant_of` | The species it is a variant of. It must have the same formula. |
| component | `processes` | Changes to the copied processes, keyed by the parent's process names. |
| process | `kind: mutation`, `from`, `to`, `per_division` | A mutation from a parent to its variant. `per_division` is a fraction, so it has no unit. |
| domain | `spreading.mechanism` | `continuum` or `cellular_automaton`. |
| domain | `spreading.interval_h` | How often biomass is spread. |

`marse check` will report the packing in each voxel at the start, the
largest growth per interval that the initial rates imply, and every species'
role with the rule it was checked against.

## Increments, and the criteria set for each

Each increment is one pull request, with its equations in theory.md, its
tests, its parameters graded in parameters.md, and a CHANGELOG entry, as the
roadmap requires of a new mechanism.

### 2d.1 Oxygen roles

This increment changes no engine code, only what a network may declare.

*Done.* What changed from the plan, each found while writing the rules:

- **Rules on factors apply to processes with a rate.** A network written for
  `marse check` alone has none.
- **An obligate anaerobe needs an inhibition factor only in a network with
  oxygen.** Without an oxygen component there is nothing for the factor to
  name, and nothing to inhibit.
- **A role is refused on a dissolved component.** It is not required on a
  particulate component that no growth process forms, such as reversibly
  bound cells or retained food.

| | Criterion | Threshold |
|---|---|---|
| D1 | Every role refuses each of its violations, and the message names the process and the factor | All 5 roles × each rule |
| D2 | An `obligate_anaerobe` grows at zero oxygen in a closed box, at the rate its other factors give | Exact, to the integrator's tolerance |
| D3 | The examples gain roles, and their results do not change | Every final digest unchanged, bit for bit |

### 2d.2 Shared space, and spreading in a column

This increment adds densities, the capacity check, the continuum provider in
1-D, the provider contract, the coupling to the solver and the manifest
entries.

| | Criterion | Threshold |
|---|---|---|
| D4 | C, N and e⁻ conserved over 10⁴ spreading steps, counting imports; a planted leak in a transfer caught at the first step | 1e-12 of the totals |
| D5 | φ ≤ 1 in every voxel after every spreading step, for every species summed (defect 8) | 1 + 1e-12 |
| D6 | The 1-D pressure flux equals the excess summed below each face (Wanner–Gujer) | To rounding |
| D7 | The splitting converges at first order in the interval. The default interval is the largest at which the thickening rate of D9 is within 0.5% of its converged value | Observed order ≥ 0.9 |
| D8 | With unlimited growth at rate μ, a neutral band starting at height z₀ is found at z₀e<sup>μt</sup>. Its error falls as the voxels shrink, since the sweep is donor-cell and smears an interface | After 3 doublings, within one voxel at 1 µm voxels; observed order ≥ 0.9 |
| D9 | A deep film fed from the liquid thickens at a constant rate, Y·J/ρ, where J is the substrate flux that the validated 1-D steady solver (V3) gives at that thickness | Constant to 2% over the last half of the run; Y·J/ρ to 1% |
| D10 | Never negative and nothing clipped, at intervals from 0.001 h to 10 h (defect 5) | Every value ≥ 0, with the ledger as in D4 |
| D11 | A run without `spreading`, and all of 2c, E1 and S1, unchanged | Every final digest unchanged, bit for bit |
| D12 | Overfilled initial states refused at load; a run reaching the top layer stopped with an error naming the time | Tests |
| D13 | Replay bit for bit, at 1 and 4 threads | Digest equal |

*Done.* What changed from the plan:

- **The sweep keeps material in order.** The well-mixed (donor-cell) sweep
  planned here failed D8: in a film growing without limit, a labelled lower
  layer reached the top of the film, and the band between the labels sat 5 to
  8 µm low without converging as the voxels shrank. In the column, what leaves
  a voxel through a face is the material nearest that face, and D8 is met.
  How this extends to faces in three directions is the first question of 2d.3.
- **The default interval is 0.25 h,** half what the D7 rule gives for the fed
  film, so that faster species keep within 0.5%. Mixing grows with the number
  of spreads, which is a second reason not to spread more often than needed.
- **D9 compares two depths, not two halves of a run.** The liquid above a
  growing film thins, so its rate rises slowly. That a deep film's rate does
  not depend on its depth is the signature the criterion was after.
- **Spreading under a salivary film is refused for now.** It arrives when the
  oral scenes need it (S2), with the film's lowest voxel as the ceiling.
- **The first `structure.csv`** holds the biovolume, the maximum thickness and
  the fullest voxel. 2d.3 adds the rest of D18's metrics.
- **One spreading block with Stage S2.** Stage S2 built its own column
  spreading for plaque at the same time: packing from the substratum up, with
  a maximum height, wear, brushing and a film riding on the plaque. The two
  were merged into one `spreading` block:
  - `mechanism: packed` is S2's (`displacement_1d_v1`), and `continuum` is
    this increment's;
  - densities live on the components, and `carried` works for both;
  - every S2 and 2d.2 result is unchanged.

  The packed mechanism settles where plaque shrinks, and the continuum leaves
  voids, so they differ only where biomass is lost. Packed spreading stays
  one-dimensional. 2d.3 extends the continuum, and 2d.5's comparison of
  mechanisms can include it in a column.

### 2d.3 Spreading in two and three dimensions

This increment adds the Poisson solve on the full region, the structural
metrics, and a 3-D example that spreads.

| | Criterion | Threshold |
|---|---|---|
| D14 | A laterally uniform 3-D box spreads exactly as the column | 1e-10 of each component's peak |
| D15 | Reflecting or swapping lateral axes moves the result with them. Reordering the components in the file changes nothing beyond rounding (defect 7, in space) | 1e-12 |
| D16 | The pressure solve agrees with a dense direct solve on small boxes, and converges to a relative residual of 1e-10 in a bounded number of iterations on irregular regions | 1e-12; iterations do not grow with the grid |
| D17 | A single colony fed from the liquid spreads with a radius that grows linearly once it is wider than three penetration depths, in a 2-D slice and in 3-D. This is the signature of [theory.md §6.2](theory.md#62-colony-radial-expansion): exponential growth of the radius would mean the coupling is wrong | Linear fit R² ≥ 0.999 over the second half; an exponential fit worse |
| D18 | Biovolume, mean and maximum thickness, substratum coverage, roughness coefficient and surface-to-volume ratio, named as COMSTAT names them ([modeling-landscape.md §5.2](modeling-landscape.md#52-structural-metrics)), computed exactly on a slab, a hemisphere and a checkerboard; written as `structure.csv` | Exact on the test shapes |
| D19 | Spreading's share of the wall time on the 3-D example, measured by profile | ≤ 20% |

### 2d.4 Variants and mutation

| | Criterion | Threshold |
|---|---|---|
| D20 | A neutral variant and its parent together evolve exactly as the parent alone | 1e-12 relative |
| D21 | Mutation in a closed box matches the closed form above | 1e-6 relative |
| D22 | A variant is exactly zero, not just small, in every voxel its cells never reached, and mutation is zero where its parent is absent (defect 6) | 0.0 |
| D23 | Each variant's share of each species, over time, written to `lineages.csv` and summarised in the manifest | Tests |

### 2d.5 The cellular automaton, and the sensitivity check

| | Criterion | Threshold |
|---|---|---|
| D24 | The cellular automaton meets D4, D5, D10, D13 and the symmetries of D15. Its tie-breaks come from the seed, and replay reproduces them | As those criteria |
| D25 | In a column, both mechanisms give the same profile of φ | To rounding |
| D26 | Two neutral colonies side by side, under each mechanism: the mixing at their boundary, as the width over which one lineage gives way to the other, is reported for both | Reported. The literature predicts more mixing under the automaton; the result is recorded either way |
| D27 | `examples/spreading_sensitivity.py` runs a two-species scene under both mechanisms and prints the D18 metrics and each species' share side by side | Runs in under 5 minutes |

## What 2d does not do

- **Detachment, erosion and sloughing.** These are in Stage 3, with dead
  biomass and the matrix. Until then a growing film needs a box tall enough
  to hold it.
- **A diffusivity that depends on biomass.** Solutes diffuse at the same
  rate in a full voxel as in liquid. Porosity and effective diffusivity are
  Stage 3.
- **Motility and chemotaxis.** A swimming cell moves relative to the biomass
  around it, which the single velocity here does not allow. They come after
  the core, as their own provider.
- **Individual cells, and random first mutants**, as decided above.
- **The rebuilt examples and periodontal study.** That is 2e, which also
  retires the version 1 expected-failure tests.
- **The performance target 2c missed.** That is Stage 7. 2d only promises
  not to make it worse (D19).

## Risks and open questions

- **Multigrid on an irregular region.** The existing V-cycle coarsens a full
  box. On Ω, the Poisson problem is symmetric positive definite, so the plan
  is conjugate gradients, preconditioned by the existing cycle on the whole
  box with p held at zero outside Ω. If D16 fails, the fallback is GMRES,
  which `marse.spatial.multigrid` already has.
- **Numerical mixing.** The donor-cell sweep smears interfaces, and mixing
  is exactly what modeling-landscape warns changes conclusions. D8 and D26
  measure it. If it dominates, a limited second-order reconstruction can be
  added to the sweep without losing positivity, as a new provider version.
- **The interval.** First order in the interval may force short intervals in
  thick films. D7 decides. An adaptive interval, bounded by the largest
  growth per voxel, is the fallback.
- **Densities.** No graded value exists yet. BM1's density is the first
  candidate. The examples' 800 C-mol per m³ colonies would overfill a
  300 C-mol per m³ voxel, and 2e has to resolve that.
- **The facultative rule may be too loose.** It does not require oxygen to
  switch the anaerobic path off. Lactic streptococci ferment in air, so
  requiring an inhibition factor would misdescribe them. This should be
  revisited with S2's five species.
- **References.** Wanner and Gujer (1986), Alpkvist and Klapper (2007) and
  BM1's density are cited from secondary knowledge. Each must be checked
  against the primary text before theory.md cites it, under the convention of
  [parameters.md](parameters.md#how-to-read-this-table).

## What changes where

- `src/marse/schemas/network.py`: densities, roles, variants and the
  mutation process.
- `src/marse/schemas/domain.py`: the `spreading` block, and the refusals of
  D12.
- `src/marse/spatial/spreading.py` (new): the providers, the Poisson solve
  and the sweep.
- `src/marse/core/reactive_transport.py`: the alternation with spreading,
  the capacity check, and `structure.csv` and `lineages.csv`.
- `src/marse/microbes/kinetics.py`: the mutation rate and its Jacobian.
- `docs/theory.md`: §6.1 is rewritten from "open question" to the model
  above, with a new section for the sweep in §9, and the assumptions in §10.
- `docs/networks.md`: the new fields, and "What is not in a network yet".
- `docs/validation.md`, `docs/parameters.md`, `docs/roadmap.md`,
  `CHANGELOG.md`.

## References

1. Alpkvist, E. & Klapper, I. (2007) A multidimensional multispecies
   continuum model for heterogeneous biofilm development. *Bulletin of
   Mathematical Biology* **69**:765–789.
2. Picioreanu, C., van Loosdrecht, M.C.M. & Heijnen, J.J. (1998)
   Mathematical modeling of biofilm structure with a hybrid
   differential–discrete cellular automaton approach. *Biotechnology and
   Bioengineering* **58**:101–116.
3. Wanner, O. & Gujer, W. (1986) A multispecies biofilm model.
   *Biotechnology and Bioengineering* **28**:314–328.
4. Wanner, O., Eberl, H., Morgenroth, E. *et al.* (2006) *Mathematical
   Modeling of Biofilms*. IWA Scientific and Technical Report No. 18.
