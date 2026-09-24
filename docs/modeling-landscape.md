# How the field models microbial growth

A survey of established practice in microbial and biofilm simulation, written
to inform the open design decisions in
[`docs/specification.md`](specification.md). It covers what the major
simulators do, how the field validates them, and — most consequentially for
MARSE — how growth in **natural environments** differs from the laboratory
growth that most kinetic parameters describe.

> **Verification.** Bibliographic records below (authors, journal, volume,
> pages, DOI) were confirmed. Numerical values and specific claims come from
> abstracts and secondary sources; publisher sites were not reachable during
> preparation, so treat these as **B**-grade under the scale in
> [`docs/parameters.md`](parameters.md) and check the primary text before
> citing any of it in a paper.

---

## 1. The established simulators

| Tool | Approach | Mechanics | Scale reported |
|---|---|---|---|
| **Wanner–Gujer** (1986) | 1-D continuum | Displacement velocity from mass balance | Analytical / 1-D |
| **Picioreanu et al.** (1998) | Hybrid differential–discrete | Cellular automaton redistribution | 2-D and 3-D |
| **BacSim / iDynoMiCS** (2001, 2011) | Individual-based | Shoving to remove overlap | Thousands of agents |
| **CellModeller** (2012) | Individual-based, GPU | Rod-shaped cells, biophysical | >30 000 cells (~100 µm colony) in ~30 min |
| **iDynoMiCS 2.0** (2024) | Individual-based, modular | **Force-based**, attractive forces, non-spherical shapes | Up to ~10 million agents in 3-D |
| **NUFEB** (2019) | Individual-based, parallel | **Discrete element method** on LAMMPS, Newtonian | Massively parallel |
| **COMETS** (2014, 2021) | Dynamic FBA on a lattice | Genome-scale metabolism + diffusion | Multi-species communities |

Two structural observations.

**The field converged on individual-based methods, then on real mechanics.**
The trajectory runs from continuum (Wanner–Gujer) through cellular automata
(Picioreanu) to individual agents (Kreft), and most recently from *shoving* —
a geometric rule that displaces overlapping spheres — to genuine **force-based
mechanics**. iDynoMiCS 2.0 made this move explicitly, and it is what enables
attractive forces and non-spherical cells, neither of which the shoving
algorithm can express. NUFEB reached the same place from a different
direction, by building on a molecular-dynamics engine (LAMMPS) and treating
mechanical relaxation as a discrete-element problem: when growth and division
push the system out of mechanical equilibrium, Newtonian equations of motion
relax it back.

**Scale is no longer the binding constraint.** Ten million agents in 3-D, or
GPU-accelerated colonies, mean the limiting factor in this field is now
parameter quality and validation, not compute. That is worth internalising
before optimising anything.

---

## 2. The biomass spreading decision

This is MARSE's main open Phase 0 question, and the literature answers it more
sharply than I expected.

Comparative work with iDynoMiCS 2.0 ranks the mechanisms by how much mixing
they produce:

$$
\text{force-based} \;<\; \text{shoving} \;<\; \text{cellular automata}
$$

and reports that **even subtle differences between spreading mechanisms
substantially affect biofilm pattern formation and population fitness** — with
the explicit recommendation that this deserves wider recognition.

This is a strong finding with an uncomfortable implication: *biomass spreading
is not an implementation detail.* It is a modelling assumption that changes
scientific conclusions, particularly about competition and cooperation, where
the outcome depends on how far related cells get separated from each other.
A cellular automaton mixes lineages more than force-based mechanics does, so
the same kinetic parameters can yield different answers about whether
cooperation is favoured.

**What MARSE should do.** Not pick one. Treat the spreading mechanism as a
declared, swappable provider — exactly the provider pattern already in
[`docs/architecture.md`](architecture.md) — and require it in the run manifest.
Any published MARSE result that depends on spatial structure should report a
sensitivity check across at least two mechanisms. This is more work than
choosing a winner, and it is the honest option given the evidence.

---

## 3. Natural states versus laboratory states

The largest gap between MARSE's current parameter set and natural conditions.
Nearly every rate constant in [`docs/parameters.md`](parameters.md) was
measured on fast-growing cells in rich or defined media. Microbes in nature
mostly do not live that way.

### 3.1 Growth rates span nine orders of magnitude

| Setting | Doubling time |
|---|---|
| Laboratory, rich medium | minutes to tens of minutes |
| Oligotrophic marine | days |
| Deep subsurface | years |

*B. subtilis* in an **oligotrophic growth state** slows from roughly 40 minutes
to roughly **4 days** — more than a hundredfold — and the nutrient level that
sustains this is comparable to *ten-thousand-fold diluted* rich medium. This
is not a stressed, dying population; it is a distinct physiological state.

There is also a sampling problem: analysis of codon usage across more than
200 000 organisms indicates **culture collections are strongly biased toward
fast growers**. The organisms we have kinetic parameters for are, by
construction, the unrepresentative ones.

### 3.2 The half-saturation constant is not a property of a species

The most actionable single finding. For *E. coli* on glucose, long-term
physiological adaptation changes $K_S$ from about **5 mg L⁻¹ to about
30 µg L⁻¹** — a factor of roughly 170, in the same organism on the same
substrate.

MARSE's [`parameters.md`](parameters.md) already warns that reported $K_S$
values span orders of magnitude. This sharpens the warning: the spread is not
only between laboratories, it is **within a single strain depending on its
history**. Treating $K_S$ as a fixed input is a modelling error, not merely an
uncertainty.

### 3.3 Cells in nature eat many things at once

Egli's work on life at low substrate concentration describes a **"multivorous"**
strategy: simultaneous uptake and metabolism of dozens of carbon substrates.
The consequences are quantitative, not just qualitative:

- Mixed-substrate growth gives **lower residual substrate concentrations** than
  single-substrate growth.
- Combining substrates **lowers the threshold concentration** at which a
  substrate can still be used.
- Consequently **generalists competitively exclude specialists** under
  mixed-substrate limitation, because the specialist cannot sustain itself at
  the residual level the generalist maintains.

Single-substrate Monod kinetics — the standard, and what MARSE implements —
cannot represent any of this. It is the right model for a controlled
laboratory experiment and the wrong one for a natural community.

### 3.4 There is a floor: the minimum substrate concentration

Below a **minimum substrate concentration** $S_{\min}$, substrate flux only
meets maintenance demand and net growth is zero. This matters "especially
during slow growth at low substrate concentrations as they are normally found
in the environment", and it has been shown to decide competition between an
enteric organism and an environmentally abundant one in glucose-limited
continuous culture.

**This is a concrete defect in MARSE's current model.** As
[`theory.md` §1.3](theory.md#13-monod-kinetics) already notes, Monod gives
$\mu > 0$ for any $S > 0$ — it approaches zero linearly but never reaches it.
In a natural-conditions simulation that is wrong: it predicts slow growth
where reality gives none, and a long simulation will accumulate biomass that
should not exist. The literature supplies both the name for the missing
concept and the mechanism: maintenance (Pirt, already implemented) sets the
floor.

A Pirt-consistent form makes the threshold explicit:

$$
\mu_{\text{net}}(S) = \mu_{\max}\frac{S}{K_S + S} - b,
\qquad
S_{\min} = \frac{K_S\, b}{\mu_{\max} - b},
$$

where $b$ is the maintenance-plus-decay rate. Note this is algebraically the
same expression as the chemostat break-even concentration in
[`theory.md` §7.1](theory.md#71-the-chemostat-and-break-even-concentration),
with maintenance in place of dilution — the two are the same idea, and MARSE
already has the code for it.

A further caution: chemostats **cannot** reproduce near-zero growth, because
growth must balance washout. Retentostats can. So parameters describing the
near-zero-growth regime that matters in nature cannot come from the standard
chemostat literature at all.

### 3.5 Dormancy and seed banks

Dormancy is a **bet-hedging strategy** producing a seed bank of individuals
that can be resuscitated when conditions improve (Lennon & Jones 2011). It is
not a modelling nicety: in soil biogeochemistry, adding explicit **active and
dormant pools** to the MEND model measurably improved validation against
270-day incubations with isotopically labelled substrates. The governing
parameters are the maximum specific growth and maintenance rates of active
microbes, and crucially **the ratio of dormant to active maintenance rates**.

This maps directly onto MARSE's adaptation layer
([`theory.md` §8](theory.md#8-adaptation-as-explicit-state-transitions)), which
already represents phenotype as transitions between declared states. Dormancy
is a well-evidenced state to include, with a literature precedent for the
parameters it needs. It also matters for a reason already noted in theory.md:
a dormant population and a population growing and dying at equal rates look
identical in aggregate biomass and behave completely differently under
perturbation.

---

## 4. Two alternative frameworks worth knowing

**Dynamic Energy Budget (DEB) theory** splits an organism into *reserve* and
*structure*. Reserves require no maintenance and provide **metabolic memory**,
smoothing fluctuations in substrate availability and allowing chemical
composition to change with conditions. A fixed fraction $\kappa$ of mobilised
energy goes to maintenance and growth. DEB is explicitly compatible with
thermodynamics, which is what lets it resolve the growth-rate/growth-yield
trade-off that Monod-plus-fixed-yield cannot.

For MARSE this is the principled answer to a problem already documented in
[`theory.md` §3.2](theory.md#32-maintenance-pirt): observed yield falls as
growth slows. Pirt maintenance handles this empirically; DEB handles it
structurally. DEB is heavier and harder to parameterise, so it belongs on the
provider-interface roadmap rather than in v1.0.

**Dynamic flux balance analysis (COMETS)** couples genome-scale metabolic
models to diffusion on a lattice. It predicted, and then experimentally
confirmed, the equilibrium species ratio of a two-member mutualistic consortium
and of a three-member engineered community. This is the strongest available
demonstration that mechanistic community models can make quantitative,
falsifiable predictions — and it is what MARSE's deferred "genome-scale
metabolism" extension point should eventually connect to.

---

## 5. How the field validates

MARSE's validation plan should adopt existing standards rather than invent
private ones.

### 5.1 The IWA benchmark problems

The IWA Task Group on Biofilm Modeling defined three benchmarks specifically to
compare modelling approaches:

| Benchmark | Tests |
|---|---|
| **BM1** | Monospecies biofilm in a completely mixed reactor: substrate flux and concentration at fixed biomass and density |
| **BM2** | Convective transport influenced by the hydrodynamics of the liquid in contact with the biofilm |
| **BM3** | Multispecies, multisubstrate: heterotrophs, nitrifiers and inert biomass coexisting |

Nine different one-dimensional models submitted solutions to BM3. These are
**published problems with published reference solutions**, which is exactly
what an external reviewer will look for.

**Recommendation.** Add BM1 and BM3 to MARSE's validation suite as V9 and V10
once the Phase 1 kernel exists. BM1 is close to the existing V3; BM3 is a
genuinely multispecies test MARSE currently has no equivalent of. BM2 requires
hydrodynamics, which is outside the v1.0 scope.

### 5.2 Structural metrics

For comparing simulated structure with confocal microscopy, the field has a
settled vocabulary, established by COMSTAT and used since: **biovolume, mean
and maximum thickness, substratum coverage, roughness coefficient, and
volume-to-surface ratio**.

MARSE should emit exactly these, under these names. Validation case V4
(attachment and biofilm initiation) currently has no quantitative target;
roughness coefficient and substratum coverage give it one that is directly
comparable with published imaging data.

---

## 6. What this changes for MARSE

Concrete consequences, in priority order:

1. **Spreading mechanism is a provider, not a choice.** Support at least
   force-based and cellular-automaton spreading; require the mechanism in the
   manifest; expect sensitivity analysis across mechanisms for any spatial
   conclusion. *(Resolves the main Phase 0 open question — by declining to
   resolve it, for good reason.)*
2. **Add a growth threshold.** Implement $S_{\min}$ via net growth with
   maintenance, so that populations stop growing where real ones do. The
   algebra is already present as the chemostat break-even calculation.
3. **Add dormancy as a phenotype state**, with the dormant-to-active
   maintenance ratio as its key parameter, following MEND's precedent.
4. **Mark every kinetic parameter with its physiological context.** $K_S$
   varies ~170-fold with adaptation history in one strain; `parameters.md`
   should record adaptation state, not just medium and temperature.
5. **Support multiplicative multi-substrate uptake** as a step toward
   mixed-substrate reality, and state plainly that single-substrate Monod does
   not describe oligotrophic communities.
6. **Adopt IWA BM1 and BM3** as validation cases, and **COMSTAT metrics** as
   the structural output vocabulary.
7. **Keep DEB and dynamic FBA as documented extension points**, not v1.0 work.

Items 2, 4 and 6 are cheap and can land before the Phase 1 kernel. Item 1 is a
design decision that shapes Phase 1. Item 3 belongs to Phase 4.

A closing note on scope. Most of this section pushes MARSE toward natural
conditions, and v1.0 is deliberately scoped to defined, laboratory-like
conditions ([`specification.md`](specification.md)). That scope remains
correct — it is what can be validated. The point of this survey is that the
*limitations* section must say clearly which of these natural-state phenomena
v1.0 does not represent, so nobody mistakes a well-mixed, single-substrate,
fast-growth simulation for a model of a soil or a marine community.

---

## References

1. Egli, T. (2010) How to live at very low substrate concentration. *Water Research* **44**:4826–4837. [doi:10.1016/j.watres.2010.07.023](https://doi.org/10.1016/j.watres.2010.07.023)
2. Harcombe, W.R., Riehl, W.J., Dukovski, I. *et al.* (2014) Metabolic resource allocation in individual microbes determines ecosystem interactions and spatial dynamics. *Cell Reports* **7**:1104–1115. [doi:10.1016/j.celrep.2014.03.070](https://doi.org/10.1016/j.celrep.2014.03.070)
3. Heydorn, A., Nielsen, A.T., Hentzer, M. *et al.* (2000) Quantification of biofilm structures by the novel computer program COMSTAT. *Microbiology* **146**:2395–2407. [doi:10.1099/00221287-146-10-2395](https://doi.org/10.1099/00221287-146-10-2395)
4. Kovárová-Kovar, K. & Egli, T. (1998) Growth kinetics of suspended microbial cells: from single-substrate-controlled growth to mixed-substrate kinetics. *Microbiology and Molecular Biology Reviews* **62**:646–666. [doi:10.1128/mmbr.62.3.646-666.1998](https://doi.org/10.1128/mmbr.62.3.646-666.1998)
5. Kreft, J.-U., Picioreanu, C., Wimpenny, J.W.T. & van Loosdrecht, M.C.M. (2001) Individual-based modelling of biofilms. *Microbiology* **147**:2897–2912. [doi:10.1099/00221287-147-11-2897](https://doi.org/10.1099/00221287-147-11-2897)
6. Lardon, L.A., Merkey, B.V., Martins, S. *et al.* (2011) iDynoMiCS: next-generation individual-based modelling of biofilms. *Environmental Microbiology* **13**:2416–2434. [doi:10.1111/j.1462-2920.2011.02414.x](https://doi.org/10.1111/j.1462-2920.2011.02414.x)
7. Lennon, J.T. & Jones, S.E. (2011) Microbial seed banks: the ecological and evolutionary implications of dormancy. *Nature Reviews Microbiology* **9**:119–130. [doi:10.1038/nrmicro2504](https://doi.org/10.1038/nrmicro2504)
8. Li, B., Taniguchi, D., Gedara, J.P. *et al.* (2019) NUFEB: a massively parallel simulator for individual-based modelling of microbial communities. *PLOS Computational Biology* **15**:e1007125. [doi:10.1371/journal.pcbi.1007125](https://doi.org/10.1371/journal.pcbi.1007125)
9. Picioreanu, C., van Loosdrecht, M.C.M. & Heijnen, J.J. (1998) Mathematical modeling of biofilm structure with a hybrid differential-discrete cellular automaton approach. *Biotechnology and Bioengineering* **58**:101–116. [doi:10.1002/(SICI)1097-0290(19980405)58:1<101::AID-BIT11>3.0.CO;2-M](https://doi.org/10.1002/(SICI)1097-0290(19980405)58:1%3C101::AID-BIT11%3E3.0.CO;2-M)
10. Rudge, T.J., Steiner, P.J., Phillips, A. & Haseloff, J. (2012) Computational modeling of synthetic microbial biofilms. *ACS Synthetic Biology* **1**:345–352. [doi:10.1021/sb300031n](https://doi.org/10.1021/sb300031n)
11. Wang, G., Jagadamma, S., Mayes, M.A. *et al.* (2015) Microbial dormancy improves development and experimental validation of ecosystem model. *The ISME Journal* **9**:226–237. [doi:10.1038/ismej.2014.120](https://doi.org/10.1038/ismej.2014.120)
12. Wanner, O., Eberl, H., Morgenroth, E. *et al.* (2006) *Mathematical Modeling of Biofilms*. IWA Task Group on Biofilm Modeling, Scientific and Technical Report No. 18. IWA Publishing.
13. Weissman, J.L., Hou, S. & Fuhrman, J.A. (2021) Estimating maximal microbial growth rates from cultures, metagenomes, and single cells via codon usage patterns. *PNAS* **118**:e2016810118. [doi:10.1073/pnas.2016810118](https://doi.org/10.1073/pnas.2016810118)

Benchmark-problem results (BM1–BM3) appear in *Water Science & Technology*
**49**(11–12), 2004, the IWA Task Group issue. Comparative results for
biomass-spreading mechanisms appear in the iDynoMiCS 2.0 paper, *PLOS
Computational Biology* (2024), [doi:10.1371/journal.pcbi.1011303](https://doi.org/10.1371/journal.pcbi.1011303).
