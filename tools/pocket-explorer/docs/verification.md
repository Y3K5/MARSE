# Package verification

Checked on an Apple Silicon Mac with Blender 5.2.2 LTS and the pinned Node
package dependencies. Geometry is uncalibrated and expert anatomical review
is **pending**. This records software and display checks, not biological or
clinical validation. The sections below record the original package build,
whose ragged basal gingival edge failed visual acceptance, and the later art
pass that replaced it (last section).

## Source package rebuild

- Rebuilt in a separate package without the bulk scaffold or full historical trajectory.
- Type checking passed; **20 viewer/geometry tests passed**, none failed.
- All 28 generated parts had zero nonmanifold edges; dentin and pulp are connected. Sampled sections and quarter cuts closed.
- Packed vertices/indices match each newly exported mesh; a cube section retains exactly three quarters of its volume.
- Preserved endpoint file/array hashes and source provenance matched. All 41,472 original endpoint values were compared before packaging; the package keeps the independently recorded extracted identities.
- Blender Boolean triangulation differed slightly from the originating export in six supporting parts. The build does not promise identical mesh bytes; exact microbial fields are unaffected.

## Exact branch browser build

The generated branch bundle rendered via local HTTP in WebGL2 at 1280 × 720,
with no horizontal overflow, zero open quarter-section chains and no observed
console warnings/errors. All nine saved taxon image identities matched the
original viewer. Healthy returned “Location unresolved” for the selected
unsupported taxon and retained a closed shallow section. The source-package
build therefore preserves the saved grid rather than projecting it onto anatomy.

The originating local art pass separately checked 390/768/1440-pixel layouts,
13 tissue selectors, keyboard orbit/Home reset, camera return and forced 2D
fallback. Those broader checks belong to that snapshot; they were not all
repeated on the package build. Earlier display execution measured about 60
frames/s on the same Mac; sustained/GPU profiling is not established.

## MARSE repository checks

- Normal Python suite: **1,037 passed**, 21 slow tests deselected, 11 existing expected failures.
- Separate slow suite: **21 passed**, 1,048 non-slow cases deselected.
- Pinned pre-commit checks passed, including ruff, gitleaks, repository structure, privacy files and Git identity.
- No engine source, protocol or simulation trajectory was changed. The package adds no Python runtime dependencies and no feedback into solver outputs.

The 11 expected failures describe existing ecosystem limitations. Passing this
suite does not turn the teaching fields into calibrated periodontal physiology.

## Publication boundary

Only 33 explicitly selected source, small synthetic fixture, screenshot and
Markdown files are proposed for Git. Bulk meshes, generated bundles, Blender
files, path-bearing Cycles metadata and logs remain outside the tracked set.
Files are below MARSE's 1 MiB size limit. Generic privacy/secret scans passed;
no private denylist is configured, so those checks cannot establish that every
possible sensitive string is absent. Reviewed screenshots have no EXIF/XMP.

The current build is an educational source contribution, not a native-app
update or a release. Rebuild requires Node dependencies and installed Blender;
disconnected installation, direct file-URL startup and other operating systems
remain unverified. No main-branch merge or deployment is included.

## Art pass on claude/pocket-art-pass (cloud verification)

Checked in a Linux cloud container, not on a Mac: 4 virtual CPUs (Intel Xeon
at 2.30 GHz under KVM), 15 GiB memory, no GPU. Browser checks used headless
Chromium 141 through Playwright 1.56.1, with WebGL2 provided by ANGLE on the
SwiftShader **software** rasterizer. Node 22.22.0 and the pinned package
dependencies. The Blender download host was blocked by the container's network
policy, so the documented authoring script ran under the official `bpy==5.2.2`
module (Blender 5.2.2 LTS, Python 3.13), imported before the script runs. This
records software and display checks, not biological or clinical validation.
Expert anatomical review is **pending**.

### Geometry rebuild

- The Blender script and Node build completed; all parts of both presets are
  closed. Every **exported triangle** edge is shared by exactly two triangles.
  The starting asset had 3–15 such defects in the ligament, trabecular bone,
  dentin and gum, from loop triangulation of Boolean n-gons. Faces whose
  triangulation would duplicate an edge are now fanned from a centre point.
- Basal gingival edge regularity (lowest gingival vertex per 2° of azimuth,
  maximum local second difference): **1.151 → 0.012** Healthy and
  **0.804 → 0.128** Periodontitis (model units). A new test fails on the starting
  asset and passes on this one.
- Soft tissue is no longer voxel-remeshed, and the cortical plates are no longer
  remeshed. Detached Boolean fragments removed from gum/connective tissue total
  at most 1.2 × 10⁻⁵ of the tissue volume; the build records them under
  `authoring_repairs`.
- Bone-top authoring changed: each azimuth keeps its crest height toward the
  tooth, and the interradicular height is −3.3 (Healthy) / −4.2 (Periodontitis).
  The earlier single central point at −7.25 made periodontitis bone slope down
  toward the tooth on every side. The selected-site margin, attachment and crest
  presets are unchanged.

### Package tests

- Type checking passed; **23 tests passed**, none failed. That is the 20 earlier
  tests plus three new ones: basal-edge regularity, exported-triangle closure,
  and plaque compartments on both cheek and tongue sides with no species
  positions.

### Browser checks (generated branch bundle, local HTTP)

- Console: no warnings or errors on any checked view. The starting build logged
  168 label-leader `Infinity` errors across the same script; non-finite
  projections are now skipped.
- Layout at 390 × 844, 768 × 1024 and 1440 × 1000: no horizontal overflow and
  no off-screen labels; WebGL2 active.
- Saved fields: in Saved model mode, all nine taxa plus "All" in both cases
  produced grid images byte-identical to the starting build (20 of 20), each
  labelled "Assumed model coordinates". Evidence tags were identical to the
  starting build: every Healthy taxon and S. sanguinis, V. parvula and
  T. denticola in Periodontitis read "Location unresolved".
- Tissue selection: all 16 selector entries and all 11 Pocket labels opened
  the matching inspector; a canvas click selected the tissue under the pointer.
- Camera and keyboard: arrows orbit; Pocket → Tooth → Pocket restored the
  remembered camera exactly; Home returned exactly to the default view; Tab
  reaches the canvas.
- State switching: Healthy/Periodontitis, cutaway/assembled, Mouth and biofilm
  detail all switched; both presets kept zero open quarter-section contours.
- Reduced motion (`prefers-reduced-motion: reduce`) loaded normally with
  `scroll-behavior: auto`; the viewer has no autonomous animation.
- Forced 2D fallback (`?render=2d`) showed the fallback drawing, disabled Mouth,
  and kept the species and saved-model controls working.
- The built-in render check measured 0.7 frames/s (starting build 0.9) at about
  545,000 triangles. That reflects CPU software rasterization only; **no GPU or
  Mac performance was measured**.

### Not repeated here

Native Blender application runs, Apple Silicon/macOS behaviour, physical
devices and real GPUs remain unverified for this pass.

## Segment pass on claude/pocket-tissue-continuity (cloud verification)

Same environment as the art pass above: a Linux cloud container with 4 vCPU,
15 GiB memory and no GPU. Headless Chromium 141 ran on the SwiftShader
**software** rasterizer, and the authoring script ran under `bpy==5.2.2`. No Mac
or GPU checks were made.

- Blender rebuild completed. Every part of both presets has 0 non-manifold edges,
  and every exported triangle edge is shared by exactly two triangles. Detached
  Boolean fragments removed from the continuous gingiva total at most 1.4 × 10⁻⁴
  of its volume; the build records them under `authoring_repairs`.
- All 768 authored 36-collar sections (192 azimuths, gingiva and core, both
  presets) are simple polygons.
- Package: type checking passed, and **24 tests passed**. New tests check the
  three-tooth segment and recorded ends, one continuous gingiva covering both
  plates along the whole segment, a smooth basal edge on each side (80 bins),
  and interdental papillae standing above the attached gingiva beside 36. The
  earlier azimuthal edge test was retired because rays toward the specimen ends
  meet the cut face; the per-side test covers the same property.
- Exported meshes now use 32-bit indices (packed format `/2`). The cortical
  plates carved by three sockets exceed 65,535 vertices.
- Browser:
  - No console warnings or errors.
  - At 390, 768 and 1440 px: no horizontal overflow and no off-screen labels.
  - In Saved model mode, 20 of 20 grid images are byte-identical to the original
    starting build, and evidence tags are identical.
  - All 16 tissue-selector entries and 11 labels open the matching inspector.
  - Healthy/Periodontitis, cutaway, Mouth and biofilm detail all switch, with
    zero open section contours, including the specimen-end trims.
  - Pocket → Tooth → Pocket restores the camera and Home resets exactly
    (isolated check).
  - Reduced motion loads normally, and the forced 2D fallback works.
- The render check measured 0.4 frames/s at about 1.06 million triangles under
  CPU software rasterization. That is not a GPU or device measurement, but the
  scene is now roughly twice as heavy as the art pass.

## Refinement pass on codex/pocket-refinement (local Mac verification)

Starting commit: `154769ac8056324139e090d9aa7fbab49fb2c1a7`, the verified
`claude/pocket-tissue-continuity` head at the start of this pass. The source
snapshot and baseline rebuild were preserved before editing. PR #54 merged during
this pass. Its merged main commit `53d5dbed5a4a94abeb3424514ba71fb76db51090`
has exactly the same Git tree as the starting snapshot; the refinement is based
on that merged commit for review. Changes are confined
to the optional viewer and this changelog; the numerical engine, protocols and
trajectories are unchanged. Expert anatomical review remains **pending**.

### Geometry and preparation

Native Blender **5.2.2 LTS**, Node **22.11.0**, pinned npm dependencies. The
complete authoring workflow and Node build completed. All **52** exported parts
have zero non-manifold edges, and every exported triangle edge has exactly two
uses. Existing closure assertions are retained. Every `buildPocket()` part closes
at both specimen ends. Packed format is `/2` with `Uint32` indices.

| Preset | Starting Pocket cutaway | Refined cutaway | Refined assembled |
|---|---:|---:|---:|
| Healthy | 1,075,273 | 386,619 | 522,401 |
| Periodontitis | 1,061,119 | 379,439 | 512,578 |

The ≤600,000 ceiling is tested for both presets and both cutaway states.
Whole exported source totals fell from 1,787,348 / 1,757,898 to
580,214 / 569,308 triangles (Healthy / Periodontitis). These counts differ from
scene counts because the viewer trims the specimen ends and adds cut faces.

Sampling was reduced before final mesh collapse. Boolean faces use the checked
export triangulation before reduction; the basal gingival rim is excluded from
collapse. Only face-free wire remnants may be removed, with a bounded assertion
and recorded repair metadata. The final rebuild needed no such wire cleanup.
Selected enamel, plaque, lumen and epithelial walls retain their display sampling.
The existing basal-edge threshold passes unchanged.

Healthy crest reference samples match the frozen starting values exactly.
Selected Periodontitis margin and attachment remain −0.15 / −6.2; distal crest
samples at x=4.2 and 4.6 remain −8.504474 / −8.506822 model units, matching the
starting authoring function. No millimetre calibration is inferred. FDI 37 is a
four-cusp context variant with cross grooves; FDI 35 is a three-cusp variant.

### Browser and Mac measurements

Actual local Mac: **Apple M1, 8 CPU cores, 16 GiB RAM, macOS 27.0.1**. Codex
in-app Chromium reports major version **155**, WebGL2 through
`ANGLE Metal Renderer: Apple M1`. This is GPU-backed browser execution, unlike
the earlier cloud software-rasterizer tests.

- Single cold Periodontitis `geometry()` observations: starting **10,622.6 ms**,
  refined **3,459.3 ms**. Healthy first preparation later measured **5,924 ms**
  with several review tabs open. These are individual observations, not controlled
  statistical benchmarks. First-use preparation can still block the UI briefly.
- Returning to a prepared state measured **0.50 ms** inside `geometry()`.
  This excludes camera framing, drawing and end-to-end interaction latency.
- The built-in short check returned **101.5 rendered frames/s** over 60 frames
  at **1175 × 540 CSS pixels**, pixel ratio 1; median CPU submission **2.00 ms**.
  This measures the browser render loop, not screen presentation, GPU duration,
  sustained performance or biological accuracy. Other hardware is untested.
- Prepared variants and anchors are cached per preset/cutaway; a test checks
  object reuse. A longer memory/performance study has not been performed.

### Interaction, evidence and responsive checks

- All **16** tissue selections opened the correct inspector. Vessel and nerve
  labels now snap to vertices of their schematic routes. The clean final browser
  session logged no errors or warnings during the checked interactions.
- Keyboard arrows changed the camera. Pocket → Tooth → Pocket restored position
  within **8.9 × 10⁻¹⁶** model units; Home reset matched exactly. Mouth, assembled
  and cutaway, Healthy/Periodontitis, biofilm detail and return controls worked.
- All nine taxon selections plus All in both presets produced **20/20** saved
  grid PNG data URLs identical to the starting build (SHA-256 comparison).
  The original 32 × 72 grid, contrast/overlap explanations, assumed coordinates
  and unresolved-location rules remain unchanged. `pocket-input.json` SHA-256 is
  `1c3485325a4738c1aab77bc26e5921d88a6b1d73d23783bf5f395ef8e2ceee4f`.
- At requested viewports **390 × 844, 768 × 1024 and 1440 × 1000**: no horizontal
  overflow and no off-stage labels, including the expanded source panel. The
  390 px view shows three default compact labels. All tissues remain selectable.
- Forced 2D fallback worked with Healthy and the species/evidence/saved-model
  controls. Mouth was disabled. Reduced-motion policy was checked using
  `?render=2d&motion=reduce`; computed scroll behaviour was `auto`. The native
  OS reduced-motion setting was not toggled. No autonomous animation is present.
- Runtime assets remain local and bundled, with no remote runtime requests or
  solver calls. A disconnected dependency installation was not attempted.

### Checks and release boundary

`npm ci`, complete native Blender authoring, `npm run build` and
`npm run check` completed: **29 tests passed**, type checking passed.
Normal engine suite: **1,037 passed, 21 slow deselected, 11 existing expected
failures**, 256.86 seconds. The slow suite was not repeated for this viewer-only
pass; earlier slow results above belong to their stated snapshots.

Pinned pre-commit checks passed, including ruff, gitleaks, repository structure,
privacy files and noreply identity. The exact outgoing-file screen found no
credentials, personal paths or private data; no private denylist is configured,
so this is not a guarantee against every possible sensitive string. The local
handoff records the screened file hashes. Rebuildable meshes, Blender/GLB files, packed buffers,
bundles, dependencies and logs stay ignored. New browser screenshots contain no
EXIF/XMP. No repository merge, deployment, native app rebuild or automatic push
is included. First-use UI blocking, stylised tissue/crown contours, schematic
oral context, unknown physical scale and pending expert review remain limitations.
