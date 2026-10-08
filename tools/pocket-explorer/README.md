# Experimental Pocket Explorer

An optional educational viewer of one generic adult human periodontal pocket
beside FDI 36. It is separate from MARSE's solver and native applications.
Authored anatomy, published qualitative observations and saved teaching fields
remain distinct inputs. Viewer controls change presentation only.

The opening view is Periodontitis / Pocket / Species / Evidence with a quarter
cutaway and labels. Mouth → Tooth → Pocket stays within one scene. Healthy
changes authored landmarks, not a progression model. All nine saved taxon
selections remain available in their original assumed 32 × 72 grid.

## Build and view

Requires Node.js 22 or newer and Blender 5.2.2 LTS for the checked authoring build.
From this directory:

```sh
npm ci
npm run anatomy
npm run build
npm run check
python3 -m http.server 8972 --bind 127.0.0.1
```

Open `http://127.0.0.1:8972/pocket.html`. Runtime scripts/textures are bundled;
there are no remote runtime assets or solver calls. Evidence links are optional
external references. Rebuild requires installed tools and dependencies; Node
installation is separate from the Python engine. The app has a 2D fallback.

Meshes, GLB, editable Blender files, packed buffers and JavaScript bundles are
rebuildable outputs ignored by Git. No bulk trajectory, raw authoring reference,
local execution log, Cycles metadata or native app is committed. The Python
script retains the earlier scaffold's reference hash but generates every live
mesh from its own rules; it does not need that reference file. No atlas triangles
were reused. This package does not add an automatic CI release or deployment.

## What goes in → what happens → what comes out

Generic Blender rules and two authored landmark presets → crown, root, pulp,
coverings, socket, separate cortical/trabecular bone and soft-tissue volumes →
an editable local anatomy and rotatable quarter cutaway. Layer widths are
exaggerated for reading and have no measured millimetre calibration.

Source-scoped human observations → schematic across-plaque detail → selected
observations and unresolved locations. Findings from selected advanced-disease
specimens never populate Healthy. Group-level observations remain groups.

Two preserved synthetic teaching endpoints → original grid and taxon controls →
read-only density images with contrast/overlap explanations. Endpoint hashes
retain the original trajectory identity as provenance; the full historical
trajectory is not included. Viewer precision is not a solver final-state digest.

## Verification and remaining limitations

The originating local art pass passed type checking and 56 software checks,
including original engine/protocol preservation checks. Those 56 are not this
package's test count. Package-specific replay results are in docs/verification.md.
All 41,472 endpoint values were compared with the original arrays before
packaging; this package locks the extracted endpoint identities independently.

**Visual acceptance remains incomplete:** the basal gingival edge over the
cortical exterior is ragged. Five selected cusps, fissures, curved taper roots,
domed furcation and connected horns/canals are authored features. Some horns
and canals lie outside a given slice. Pocket crops apical context intentionally;
Tooth shows the full model. Whole-mouth context remains earlier procedural art.
Expert anatomical review is **pending** and physical calibration is unresolved.
Blender Boolean triangulation can vary slightly across rebuilds: this package
does not promise byte-identical mesh exports. Saved endpoint data stays exact;
viewer buffers are checked against each newly exported source.

![Quarter cutaway](docs/images/cutaway.jpg)

![Exterior with remaining gum-edge defect](docs/images/exterior.jpg)

Software/topology checks establish their stated mesh/display properties only.
This is not clinical validation, measured pocket anatomy, diagnosis, or a model
of microbial contacts, oxygen gradients, immune function or disease progression.
Canine inference and immune binding/epitope scoring are excluded. MARSE's v1.0
scope is unchanged. Nothing feeds back into solver trajectories.

## Sources and licences

- [Oral Biofilm Architecture on Natural Teeth](https://pmc.ncbi.nlm.nih.gov/articles/PMC2827546/): selected human specimens and across-biofilm observations; not universal bands down a pocket.
- [Mandibular first-molar micro-CT research](https://pmc.ncbi.nlm.nih.gov/articles/PMC3662034/): qualitative root/canal morphology; no parameter calibration here.
- [Root/canal variation](https://pubmed.ncbi.nlm.nih.gov/30147408/): the chosen arrangement is not universal.
- [Human dentoalveolar micro-CT](https://pubmed.ncbi.nlm.nih.gov/29732712/): separate cortical/trabecular/socket context; not model calibration.

Viewer source and authored mesh rules are a MARSE contribution under the
repository's Apache-2.0 licence. Three.js is MIT; dependency licences remain in
the installed packages and build output. No Z-Anatomy mesh is distributed here.
