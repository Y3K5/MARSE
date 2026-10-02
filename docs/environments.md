# Environments

A microbe lives somewhere: on a surface, under a liquid, in a flow, among other
cells, near or far from food. MARSE builds these settings as **scenes**:
experiment files that say what the surface is made of, what flows over it,
which cells the liquid carries, and how those cells bind.

This page describes:

- how a scene is assembled;
- the scenes that exist, with their results;
- the programme that will add more.

The model behind binding is in [theory.md §6.4](theory.md#64-attachment-to-surfaces),
the fields in [networks.md](networks.md#cells-binding-to-the-surface), and the
parameters and their grades in [parameters.md §7](parameters.md#7-surfaces-and-adhesion).

## The programme

The environments come in increments, each planned before it is built:

| Increment | Delivers | Status |
|---|---|---|
| **E1** | Binding to flat surfaces, from the delivery of cells by the flow to reversible, then locked, binding limited by the free area. The locked cells grow. A laboratory flow chamber and a dental scene with four materials. | this release |
| E2 | Binding from surface physics, the extended DLVO theory, for surfaces and waters that have no adhesion measurements. It covers van der Waals attraction, the electrical double layer as ionic strength and pH change, and acid–base interactions. Natural waters by region (river, lake, groundwater, seawater) and rock and mineral surfaces (quartz and granite, calcite and limestone, basalt, pyrite). | planned |
| E3 | Soil as a full 3-D pore structure. Solid grains and pores, each grain of one material (quartz, clays, iron oxides, organic coatings), binding on every grain face, and water flowing through the pores. Soil types by composition. | planned |
| E4 | Geographic scenes, from measured regional compositions of soils and waters, built with E2 and E3. | planned |

A substratum face is the unit that a material, and binding with it, attaches
to. E3 gives every grain in a soil faces of its own, in every direction, so it
reuses the machinery of E1 rather than replacing it.

After E1, the roadmap turns to saliva, diet and caries. That stage compares
high- and low-sugar diets over months to years, with salivary flow, clearance
and buffering, and the pH drops that select acid-tolerant species. Its first
increment, S1, brings the mouth, the diet and the Stephan curve:
[the oral scenes](#the-oral-scenes) below. Host biology (immune, structural
and functional cells, and their signalling) comes after v1.0, as a separate,
clearly labelled layer ([roadmap.md](roadmap.md)).

## A scene

A scene is an experiment file with a domain in space
([networks.md](networks.md#running-a-network-in-space)). Five fields say
what the surface is and what binds to it, and they come together:

| Part | Field | What it states |
|---|---|---|
| The surface | `substratum` | Patches of material covering the bottom of the box, and the conditioning film that coats them. |
| The liquid | `liquid` | Its temperature and viscosity. |
| The flow | `flow` | Its shear rate at the surface, and how far downstream the scene sits. |
| The cells in the liquid | `suspension` | Each species that binds: cells per mL, cell diameter, carbon per cell, and the area a bound cell blocks. |
| How they bind | `adhesion` | For every species on every material: the attachment efficiency, and the rates of detachment and locking. |

The rest of a scene is a network in space, as in Stage 2c:

- **The resources** are the dissolved components, with their concentrations in
  the bulk liquid and their diffusivities.
- **The organisms** are the particulate components and the processes by which
  they grow.

Each binding species has two components with the same formula:

- one holds its reversibly bound cells;
- the other is its biomass, which locked cells join.

Only the biomass may take part in processes.

**Check.** `marse check` previews a scene before it runs. It prints:

- the share of each material;
- how fast each species is delivered to the surface;
- how fast it binds to each material at first;
- how long a surface would take to get half-way to jamming, if nothing
  detached or grew.

**Run.** `marse run` writes `surface.csv` beside the usual outputs. For each
material it records the bound cells per cm² of each species, reversible and
locked together, and the fraction of the material covered by cells.

"Covered" is the area under cells, which is what microscopy measures. It is
the area covered by discs of the cells' diameters placed at random,
$1 - \exp(-\sum_s n_s \pi d_s^2/4)$. The blocked area, which limits binding,
is larger: a bound cell keeps others from binding in a neighbourhood around it.

## The laboratory scene

[`examples/environments/lab/flow_chamber.json`](../examples/environments/lab/flow_chamber.json)
is a parallel-plate flow chamber, the standard instrument for measuring how
bacteria stick. *Streptococcus oralis* is suspended in buffer at 3 × 10⁸
cells per mL and flows over glass at 37 °C, at a wall shear rate of 15 s⁻¹,
20 mm from the inlet. Half the glass is bare, and half carries a salivary
conditioning film. The buffer holds no nutrients, so nothing grows, and the run
measures how many cells stick, and how fast, over four hours.

```bash
marse check examples/environments/lab/flow_chamber.json
marse run examples/environments/lab/flow_chamber.json
```

The cells diffuse at 0.73 µm² s⁻¹ and arrive at 1.19 × 10³ cm⁻² s⁻¹.
Binding then gives, per cm²:

| Glass | Binds at first | 30 min | 1 h | 2 h | 4 h | Covered at 4 h |
|---|---|---|---|---|---|---|
| Bare | 595 s⁻¹ | 7.5 × 10⁵ | 1.39 × 10⁶ | 2.43 × 10⁶ | 3.78 × 10⁶ | 2.4% |
| Saliva-coated | 357 s⁻¹ | 5.5 × 10⁵ | 1.04 × 10⁶ | 1.89 × 10⁶ | 3.12 × 10⁶ | 2.0% |

- **The initial rates** lie inside the range measured for oral streptococci on
  glass in such chambers, 0 to 2.9 × 10³ cm⁻² s⁻¹ (Sjollema, Busscher and
  Weerkamp 1988).
- **The curves bend.** Each approaches the jamming density of 5.5 × 10⁶ cm⁻²
  as bound cells block the surface, the shape flow-chamber data are fitted
  with.
- **The two halves differ** only through the attachment efficiency and
  detachment rate of each material. Those parameters are illustrative, set to
  give plausible magnitudes. Replace them with values fitted to a flow-chamber
  experiment before drawing conclusions from the scene.

## The dental scene

[`examples/environments/dental/dental_surfaces.json`](../examples/environments/dental/dental_surfaces.json)
is one surface patterned in four stripes, all under a salivary pellicle:

- enamel;
- titanium;
- zirconia;
- acrylic (PMMA).

The stripes sit side by side under the same salivary film for 24 hours. The
box is 80 × 20 × 40 µm in voxels of 5 µm.

- **The film** is 80 µm thick and moves at 2 mm per minute, within the
  measured 70–100 µm and 0.8–7.6 mm per minute. So its wall shear rate is
  $3\bar u/\delta$ = 1.25 s⁻¹.
- **The early colonizers** *S. oralis* and *S. sanguinis* are suspended in the
  saliva at 10⁷ and 5 × 10⁶ cells per mL. They bind to the pellicle and lock.
- **They grow** on salivary glucose by lactic fermentation, taking nitrogen
  from ammonium.

```bash
marse check examples/environments/dental/dental_surfaces.json
marse run examples/environments/dental/dental_surfaces.json   # 35 s
python examples/surface_adhesion.py   # the four materials as fast columns, 7 s
```

Bound cells per cm², both species together:

| Material | 1 h | 6 h | 12 h | 24 h | Covered at 24 h |
|---|---|---|---|---|---|
| Enamel | 4.9 × 10⁴ | 4.9 × 10⁵ | 1.93 × 10⁶ | 1.93 × 10⁷ | 12.0% |
| Titanium | 4.9 × 10⁴ | 4.9 × 10⁵ | 1.93 × 10⁶ | 1.93 × 10⁷ | 12.0% |
| Zirconia | 3.1 × 10⁴ | 3.1 × 10⁵ | 1.25 × 10⁶ | 1.27 × 10⁷ | 8.1% |
| Acrylic (PMMA) | 4.9 × 10⁴ | 4.9 × 10⁵ | 1.93 × 10⁶ | 1.93 × 10⁷ | 12.0% |

What the run shows, and what it does not:

- **Binding, then growth.** In the first hours the count rises by binding, at
  9.5 *S. oralis* and 4.4 *S. sanguinis* per cm² per second on enamel,
  titanium and acrylic. Later it rises by the growth of locked cells. By 24 h
  growth accounts for nearly all of it.
- **Titanium against zirconia.** Zirconia's attachment efficiency is set to
  0.63 of titanium's, calibrated on the 12.1% against 19.3% of surface covered
  after 24 h in the mouth (Scarano et al. 2004). The scene's 8.1% against
  12.0% reproduces that calibration; it is not a prediction. Another study, of
  low-roughness surfaces in the mouth, found no difference
  ([parameters.md §7.2](parameters.md#72-titanium-against-zirconia)).
- **Enamel, titanium and acrylic** share every parameter, because no
  comparison between them was found. Yet enamel ends 0.1% below the other two.
  The box is one tile of a repeating surface, so the enamel stripe lies between
  titanium and acrylic, which grow well, and competes with both for glucose.
  Titanium and acrylic each border one strong stripe and zirconia, and they
  agree.
- **Binding stops too early.** The blocked area per cell is a property of
  single cells. Cells that divide stay together in microcolonies, which block
  less area per cell than scattered cells. E1 counts every cell as blocking
  the same area. On titanium, binding therefore stops at about 17 h, when
  less than 4% of the surface is under cells. Coadhesion, cells binding to
  cells already bound, is not modelled either.
- **Glucose** is a steady daytime average. Meals, and the pH they bring, come
  with the saliva-and-diet stage.
- **Colonies grow in place.** A column can spread its biomass
  ([networks.md](networks.md#biomass-that-spreads)), but this scene is a 3-D
  box, which spreads with Stage 2d, increment 2d.3.

## Assumptions of E1

- **Surfaces are flat,** and the bottom of the box. Roughness belongs to
  retention, which a later increment models.
- **No sedimentation.** The scenes are vertical, or face down. A surface facing
  up also collects cells that settle, which E1 does not model.
- **A fixed suspension.** Binding does not deplete the cells in the bulk, and
  cells are not tracked through the liquid inside the box.
- **Stated efficiencies.** The attachment efficiency of each species on each
  material is a number in the scene, not derived from the surfaces' chemistry.
  E2 derives it.
- **Continuum cells.** Bound cells are densities, not individuals. Below about
  one cell per face the numbers are expected values.

The ecosystem model's surface transfer ([adhesion-detachment.md](adhesion-detachment.md))
is a separate, phenomenological rule of the version 1 engine, not this model.

## The oral scenes

A tooth's plaque does not sit under a well-mixed liquid, but under a thin
film of saliva that the mouth renews, and it sees sugar only when the mouth
does. The three scenes in
[`examples/environments/oral`](../examples/environments/oral) model one site:

- **The plaque** is a column of 150 µm, with one acidogenic population that
  ferments sugar to lactic acid, more slowly as the pH falls. Carboxyl groups
  of the cell walls and matrix buffer it and hold potassium, which they
  release as they take up protons ([theory.md §3.8](theory.md#38-acidbase-equilibria-and-ph)).
- **The film** is 100 µm of saliva over it, moving at 6 mm per minute and
  renewed from the mouth ([§4.8](theory.md#48-a-salivary-film-and-the-mouth)).
- **The mouth** secretes resting saliva, and more and more alkaline saliva as
  it tastes sugar, with the buffers Bardow et al. (2000) measured. It
  swallows as Dawes (1983) described.
- **The diet** is what differs ([§4.9](theory.md#49-the-diet)):
  - [`stephan_rinse.json`](../examples/environments/oral/stephan_rinse.json):
    10 mL of 10% sucrose held for a minute and spat out, the challenge of a
    Stephan curve;
  - [`sipping.json`](../examples/environments/oral/sipping.json): 100 mL of a
    drink of 10% sucrose sipped over 20 minutes;
  - [`pocket.json`](../examples/environments/oral/pocket.json): the rinse,
    after which food particles holding sugar stay on the teeth and dissolve.

```bash
marse check examples/environments/oral/stephan_rinse.json
marse run examples/environments/oral/stephan_rinse.json   # 10 s
python examples/stephan_curve.py   # the three scenes side by side, 44 s
```

`marse run` writes `ph.csv` (the pH at the substratum and its range in the
box) and `mouth.csv` (the mouth's volume, flow, swallows and composition)
beside the usual outputs. The pH at the substratum, over 90 minutes:

| Scene | Lowest pH | Minutes below pH 5.5 | Back above pH 6 |
|---|---|---|---|
| The rinse | 4.87 at 15.6 min | 31.7 | 42.9 min |
| Sipping | 4.80 at 34.8 min | 52.3 | 64.7 min |
| Food left on the teeth | 4.78 at 22.2 min | 56.2 | 69.5 min |

What the scenes show, and what they do not:

- **A Stephan curve.** After the rinse the pH falls by 1.9 units, from 6.80
  to 4.87 at 16 minutes, and is back above 6 within 45 minutes, as measured
  curves do (Stephan 1944). The acid production, the buffer, the
  plaque's thickness and the film's speed were calibrated to criteria set
  before the stage was built ([validation.md](validation.md#the-stephan-curve)),
  so this shape is a calibration, not a prediction.
- **What makes a high-sugar eater's plaque acid for longer.** Sipping the same
  sugar over 20 minutes keeps the plaque below pH 5.5 for 52 minutes instead
  of 32, and a little food left on the teeth for 56. Low salivary flow, a
  slower film, thicker plaque and more buffer move the curve the ways the
  literature reports. Those directions are predictions of the model, and the
  tests check each one.
- **One population, one site.** Which species make the acid, and how the pH
  selects among them, is Stage S2. Enamel dissolving below the critical pH is
  Stage S3, and many sites and years of diet are S4.
- **Nothing makes base.** Urea and arginine, which plaque turns into ammonia,
  are not yet modelled, so the return to neutral relies on the saliva alone.
