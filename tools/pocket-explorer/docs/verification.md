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
