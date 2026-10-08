# Package verification

Checked on an Apple Silicon Mac with Blender 5.2.2 LTS and the pinned Node
package dependencies. Geometry is uncalibrated and expert anatomical review
is **pending**. This records software and display checks, not biological or
clinical validation. Supporting-tissue visual acceptance still fails because
of the ragged basal gingival edge.

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
