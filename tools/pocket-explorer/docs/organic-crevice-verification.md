# Organic crevice verification — October 9, 2026

This pass changes authored presentation geometry, materials, normals and camera
flow in the existing Crevice scale. It remains one scene, separate from Anatomy
Studio and numerical simulation. No physical scale or measured registration is
introduced. Expert anatomical review is pending.

## Checked results

- Build and TypeScript checks pass. All 45 package tests pass. A new mesh test
  checks all nine tissue volumes in both presets: every welded edge has two
  opposing incident faces, signed volume is positive, every triangle has nonzero
  area, and exposed surfaces/section caps have separate material groups.
- The existing 20 density-display RGBA arrays are identical to the baseline.
  All 163 frozen engine/example/test/input/organic anatomy files are unchanged.
  Scientific records, morphology references and stylesheet also retain the
  immediate pre-pass hashes.
- Browser: nine selectors retain correct names, sources, nine exemplars and
  unresolved location status. Actual ray picking on focused P. gingivalis
  changes Structure to Species and opens its morphology inspector. All three
  camera views operate; manual arrow input interrupts travel and the pose
  returns exactly through Pocket → Crevice. Home resets.
- Reduced-motion query uses immediate camera changes. 2D fallback disables
  Crevice and retains Pocket. Healthy keeps location unresolved rather than
  importing disease findings. No captured warning/error logs.
- Section layouts at 390, 768 and 1440 px have no horizontal overflow or
  overlapping tissue-label rectangles. Camera controls are below the canvas.
- Apple M1 / ANGLE Metal: 60.4 frames/s over 60 frames, 1.90 ms median CPU
  submission, 569 × 570 CSS-pixel stage, pixel ratio 1.5, 130624 triangles.
  This short execution sample is not sustained profiling or biological evidence.
- Nine refreshed GLBs pass header, chunk, buffer-bound and identity/hash checks.
  This is structural verification, not complete external glTF conformance.

## Visual inspection and artifacts

The initial block silhouettes were inspected before material work. Final browser
screenshots compare Opened section, Look down and Biofilm surface at the same
1440 × 1000 viewport, Periodontitis preset, all-taxa selection and labels enabled.
Camera framing is intentionally improved rather than held at the old pose.
Healthy and mobile screenshots are also retained. A frozen baseline bundle was
compiled from preserved pre-pass source into a separate local receipt directory;
a temporary loopback-only server supplied baseline captures without reverting
live source. Early exploratory captures are not the final comparison.

The organic-crevice receipt includes comparison.html, before/after JPEGs,
package-check.log, integrity.json, browser-checks.json, object-checks.json and
microbial-objects.zip. The updated start guide explains controls and boundaries.
Everything remained local; nothing was uploaded, published, committed or pushed.

## Remaining limits

Closed seams do not prove absence of all mesh self-intersections or tissue
intersections. This is authored educational anatomy; the cut region and widths
are exaggerated. The fold/core share a geometric boundary, not measured
histology. Texture is an artistic cue, not a reconstructed fiber or trabecular
network. Cell sizes and placements are independent display choices; no natural
contacts, abundance, niches or growth behavior are inferred. Morphology alone
cannot establish species identity. Canine models and immune predictions remain
outside this pass. Clinical and anatomical expert review is pending.
