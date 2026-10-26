# Crevice anatomy refinement — October 10, 2026

This local pass rebuilds the authored Crevice geometry, its materials, labels and
three cameras. It stays one scene with one inspector, separate from Anatomy
Studio and the numerical engine. It is stylized educational anatomy with no
physical scale or measured registration. Expert anatomical review is pending.

## What changed

- **One profile swept around the root.** Every tissue is a closed outline in a
  root-to-cheek plane swept around a single tapering root axis. The gingiva
  therefore wraps the tooth as a continuous cuff. The near end is the
  deliberate section face. The apical end, the far end and the region near the
  root axis dissolve into the stage instead of ending as slabs.
- **Separately readable tissues.** Enamel thins to a feather edge at the
  cementoenamel junction. Cementum thickens apically over a tapering, slightly
  curved root with pulp and canal. Separate closed volumes cover:
  - sulcular epithelium;
  - a wedge-shaped junctional epithelium on the tooth;
  - gingival connective tissue, which meets cementum between the junctional
    epithelium and the crest;
  - oral gingival epithelium with a rounded margin;
  - periodontal ligament;
  - a rounded alveolar crest.

  Neighbouring tissues share boundary vertices. A test samples each section
  densely and confirms no point lies inside two tissues.
- **Materials.** Surfaces are mostly matte, with low sheen on soft tissue and
  no strong clearcoat. Section faces are flatter and slightly paler.
  Connective tissue, ligament and bone sections carry faint artistic pattern
  cues. The plaque film has gentle relief and mottling on its fluid side only.
- **Microbial forms.** The nine identities, their geometry, references and
  controls are unchanged (all 9 geometry and GLB hashes match the pre-pass
  build). Forms now rest on the plaque film below the retained margin, rather
  than floating in the fluid. The gallery stays labelled as illustrative.
- **Three distinct relationships:**
  - *Opened section:* layers from tooth to bone, in a three-quarter view of
    the cut face.
  - *Look down:* the cuff and its narrow entrance continuing around the root,
    with the V of the sulcus at the cut face.
  - *Biofilm surface:* plaque on the tooth. The soft tissue is cut away, while
    a margin band, the junctional attachment and the crest are retained.
- **Labels.** Labels stand in two columns outside the model silhouette and
  never on the pocket. Each leader ends on an anchor computed from the same
  profile curves as the mesh. After the camera settles, anchors hidden behind
  tissue get dashed leaders. Narrow screens use a short essential set.
- **Cameras.** Camera travel takes 0.72 s, orbits around the target and eases
  in and out. Pointer, wheel and arrow input stop it; reduced motion jumps
  immediately; Home resets the current view.

## Checked results (cloud container)

- **Typecheck and build:** both pass.
- **Package tests:** 46 of 49 pass, including four new crevice tests:
  - profiles are simple and do not overlap;
  - the epithelial, connective and ligament relationships hold;
  - section label anchors lie inside the tissue they name, on the cut face;
  - forms rest on the film below the margin band.
- **Atlas test failures:** the three failing tests are the BodyParts3D atlas
  tests. This sandbox's network policy blocks the atlas source host, so a
  clearly labelled local placeholder atlas buffer was used to compile. The
  same three tests fail on the unmodified baseline here, and Mouth/Face were
  not assessed in this pass.
- **Diff scope:** only four viewer files changed: `src/pocket-crevice.ts`,
  `src/pocket.ts`, `pocket-explorer.css` and `tests/pocket-crevice.test.mjs`.
- **Saved fields:** all 20 saved density displays are byte-identical to the
  pre-pass build, and the saved input SHA-256 is unchanged.
- **Microbial forms:** the nine microbial geometries, the claim text and the
  nine GLB exports are identical.
- **Browser checks:** 26 of 26 pass:
  - travel and settle;
  - wheel and arrow interruption, Home, and reduced motion;
  - pose restore through Pocket → Crevice;
  - label selection;
  - picking at every section anchor, which selects the named tissue (12/12);
  - leaders ending on anchors, no overlapping label boxes, and no label over
    any of 756 pocket-opening samples;
  - nine identities with sources;
  - each taxon opening its own identity with "location unresolved" and a
    focused camera;
  - Healthy keeping placements unresolved;
  - the 2D fallback disabling Crevice;
  - no overflow at 390, 768 and 1440 px;
  - no console warnings or errors.
- **Hardware:** headless Chromium 141 with the SwiftShader software renderer,
  measured in this cloud container, not a Mac:
  - baseline: 1.6–1.7 frames/s and about 3.2 ms median CPU submission over
    two runs, 130,624 triangles;
  - after: 1.6 frames/s and about 6.0–6.3 ms, 141,402 triangles.

  This shows the relative cost under software rendering only and says nothing
  about GPU performance.

## Honest remaining defects

- **Crown:** it is a smooth dome standing in for coronal context, not a molar
  crown.
- **Root shape:** the swept patch is round about one axis. Real molar roots are
  not round, and the root trunk and furcation are not represented.
- **Junctional length:** it is an authored half of the attachment-to-crest
  distance.
- **Healthy junctional position:** the healthy attachment preset sits below
  the cementoenamel junction. The junctional epithelium therefore lies on
  cementum in Healthy; presets were not changed.
- **Rete ridges:** they are an artistic wave and constant around the sweep.
  Fibre groups, lamina dura, cortical/trabecular bone and vessels are not
  modelled; bone is one volume with a section pattern.
- **Widths:** ligament width, sulcus width, plaque thickness and relief, and
  cell sizes are exaggerated or artistic. Cell size also changes between
  presets.
- **Biofilm-surface margin band:** the retained band occupies the same space
  as the full gingiva and is only shown while the full gingiva is hidden.
- **Picking:** on the gingiva, oral epithelium versus margin is chosen by a
  height threshold.
- **Labels:** occlusion dimming updates only after the camera settles. At
  390 px, label columns may overlap tissue edges, though not the pocket.
- **Mesh audit:** forms may sink into the film relief by up to about 0.045
  authoring units. There is no full 3D self-intersection audit beyond the
  closed-seam and section-overlap tests.
- **Face buffer:** it was rebuilt with the pip `bpy` 5.2.2 module, and its hash
  differs from the lock (4,256,375 vs 4,256,381 bytes). The lock was not
  edited.
- **Review:** expert anatomical review is pending. Passing software checks is
  not biological or clinical validation.

Nothing was committed, pushed, uploaded or published in this pass.
