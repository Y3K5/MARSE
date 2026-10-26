# Magnified crevice verification — October 9, 2026

This local pass adds an authored explanatory patch and nine qualitative 3D
microbial forms. It does not alter the numerical engine, saved endpoints or
organic anatomy assets. The default remains Pocket / Periodontitis / Species /
Evidence; `?view=crevice` opts into the new scale.

## Checked

- Build and TypeScript check pass. All 44 package tests pass, including four
  new tests for deterministic finite geometry, distinct forms, authored landmark
  ordering and camera poses, immutable inputs and explicit placement limits.
- All 20 density display RGBA arrays match the baseline byte for byte.
- All 163 frozen engine/example/test/input/organic anatomy files retain hashes.
- Nine GLB exports pass header, chunk, buffer-bound, identity and SHA-256 checks.
  These are structural checks, not a complete external glTF conformance audit.
- Mac browser: all nine selectors open the correct identity/source; actual 3D
  picking selected P. gingivalis; all three cameras change pose; arrow-key orbit
  returns exactly after Pocket → Crevice; Home resets. Healthy keeps positions
  unresolved and does not infer organism presence.
- Layouts at 390, 768 and 1440 pixels have no horizontal overflow. Camera
  controls remain below the rendered canvas. 2D fallback disables Crevice and
  retains Pocket and evidence/model controls. No captured warning/error logs.
- Apple M1 / ANGLE Metal: 60.4 frames/s over 60 frames, median CPU submission
  1.10 ms, 569 × 570 CSS-pixel stage, pixel ratio 1.5, 61,056 triangles.
  This short sample measures display execution, not sustained performance.

## Scientific and visual limits

Widths, cell sizes and artistic surface relief are uncalibrated. Examples are
laid out for inspection, not measured niches, contacts, abundance or succession.
Display-fit sizes differ between presets without biological meaning. Morphology
alone cannot identify a species; Actinomyces remains a group. Available source
records retain specimen context, including extraoral material where applicable.
No molecular envelope, adhesion sites, virulence functions or growth parameters
are inferred. Saved density fields remain in their own original coordinates.

The crevice is independent magnified authored geometry, not measured segmentation
or registration of the tooth. Closed profile construction is used, but this pass
has no comprehensive intersection/watertightness audit of the new tissue meshes.
This remains stylized educational anatomy. Expert review is pending; passing
software checks establishes no clinical or biological validation.

## Local receipts

The existing pocket-anatomy receipt directory contains package-check logs,
browser/layout records, protected-file comparisons, map comparisons, screenshots
and `crevice-objects.zip` with nine GLBs and their provenance manifest. No files
were uploaded, published, committed or pushed in this pass.
