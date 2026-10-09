# Local anatomy and species-map verification — 2026-10-09

This record describes software and presentation checks. It does not establish
clinical accuracy or biological predictive ability.

## Inputs and separation

The detailed pocket and its two landmark presets are authored illustrations.
Only the crown presentation uses a cited cohort size reference; it does not
calibrate surrounding tissue or pocket dimensions. Source geometry is cloned
for presentation. Original organic anatomy bytes remain unchanged.

BodyParts3D release 4.0 reduced OBJ surfaces provide the independent Face/Mouth
reference. The subset retains 63 named source elements and SHA-256 identities;
53 are used in the reviewed native scene after excluding isolated vessels and
salivary glands. The viewer preserves original float32 vertex coordinates under
a rigid basis change and an authored lower-jaw pose. There are 28 uniquely named
teeth; FDI 36 is FMA55704/FJ1254. Palatal elements are muscle surfaces, not a
complete soft-palate mucosal reconstruction. Source units are millimetres.

Morphology records contain qualitative form, source URL and specimen context.
Forms are drawn independently of density bins. Colors are identity colors;
cell number, relief, normalized size and drawn arrangement are illustrative.
The saved source is human-only, arbitrary model density, row-major 32 × 72,
with no physical anatomical registration. Evidence positions and saved grids
remain separate. No human evidence is transferred into canine parameters.

## Measured checks

- npm run check: TypeScript passed; 40 tests passed, zero failures.
- Frozen files: 163 hashes unchanged, including engine code, examples, tests,
  source saved input and original organic anatomy.
- Saved map: 20 of 20 RGBA arrays are byte-identical to baseline commit
  4abdb5df9b7e80c91bc1ecfd8f4a40886bd50932 (two presets × All/nine taxa).
- Browser: all nine taxon selections showed matching grid IDs and morphology
  source links. Healthy P. gingivalis remained “Location unresolved.”
- Keyboard: End reached tissue-response stage 5 with the bone-defect reference;
  opening the story from Structure correctly returned to Species detail.
- Forced 2D fallback retained saved fields and species forms; atlas-only camera
  buttons were disabled. Reduced-motion query loaded successfully. A sustained
  animation timing audit under reduced motion was not performed.
- Layout: no document horizontal overflow at 390, 768 and 1440 CSS pixels.
  At 768 the map stacks above the gallery rather than leaving tiny cards.
  Detail flows down the page at tablet/phone widths, avoiding nested scrolling.
- Console sample: no errors or warnings in the checked final browser sequence.
- Actual Mac Mouth sample: 101.6 frames/s over 60 rendered frames; median CPU
  submission 1.20 ms; 1040 × 755 CSS stage, pixel ratio 1, Apple M1 Metal/ANGLE,
  117,812 triangles. This short execution sample is not sustained GPU profiling.
- Runtime code/assets are bundled and source links are optional. There was no
  network-disconnected browser replay; offline packaging is a source/asset
  check, not a claim that such a replay was performed.

## Native local Blender verification

The installed local Higgsfield Blender MCP was exercised through its stdio
protocol, not a hosted scene job. It saved the editable reference scene and
rendered a 1100 × 1100 Cycles image. Six isolated-light renders were inspected.
The GLB parses as glTF 2, contains 53 meshes and 53 source-named nodes, and is
4,274,936 bytes. Procedural Blender noise/bump materials are not baked into GLB;
its basic PBR appearance differs from the native render.

Raw atlas boundaries largely represented coincident-vertex seams. A diagnostic
weld on temporary copies left three boundary edges in upper gingiva FJ1252
(one small triangular opening). Original atlas meshes were not silently repaired.
Native scene source identities and the independent pocket remain separate.

## Limits and review status

Expert anatomy review: **pending**. Pocket surfaces, visible tissue widths,
canals and nerve courses are illustrative. No calibrated pocket ruler, patient
scan, external facial skin, measured TMJ dynamics, cell-contact map, molecular
binding sites, microbial growth clock, immune protection or clinical diagnosis
is supplied. The reference image's full facial tissue fidelity is not achieved.
The full MARSE engine and slow suites were not rerun in this presentation pass;
their files were hash-checked. Nothing was committed, uploaded or pushed.
