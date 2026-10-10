# Pocket Explorer — start here

This is a local educational viewer for a human tooth–gum interface. It is
independent of Anatomy Studio and the MARSE numerical engine. Expert anatomical
review is pending.

Reference atlas surfaces, authored pocket anatomy, source-linked observations
and two saved teaching endpoints go in. The viewer changes cameras, visibility
and selected information. A rotatable mouth, readable tooth section, enlarged
density map and explanatory species forms come out. Controls do not run or
change a microbial simulation.

Open the local preview at http://127.0.0.1:8975/pocket.html while its server runs.
Choose **Crevice** for a magnified tooth–gum space, or open
http://127.0.0.1:8975/pocket.html?view=crevice directly. The crevice is one small
piece of tooth and gum, cut open at the end nearest you and curved around the
root, so the gum wraps the tooth like a cuff. Its far, deep and inner edges fade
into the background on purpose: they are where the drawing stops, not anatomy.

Each camera answers one question:

- **Opened section** — what lies between tooth and bone? Reading outward from
  the tooth, the cut face shows:
  - enamel, root dentin and pulp, with cementum on the root;
  - attached plaque and the fluid space;
  - sulcular epithelium lining the soft-tissue wall;
  - junctional epithelium attached to the tooth at the base;
  - connective tissue, which meets the root between the junctional epithelium
    and the bone;
  - periodontal ligament and the rounded alveolar crest.
- **Look down** — how does the gum meet the tooth? The gingival margin hugs the
  crown like a cuff around a narrow entrance. At the cut face the entrance opens
  into the V of the sulcus or pocket.
- **Biofilm surface** — what coats the tooth? The soft tissue is cut away to show
  plaque on the root. A margin band, the junctional attachment and the crest are
  kept for orientation. Nine selectable 3D forms rest on the plaque. Click a form
  or use the species selector to focus it and open its morphology reference.

Labels sit to either side and never cover the pocket. Click a label, or the
tissue itself, to read about it in the inspector. A dashed leader means its point
is behind tissue from the current angle. Camera buttons travel smoothly for
under a second. Dragging, scrolling or arrow keys stop the move, and Home
returns to the current view. Reduced-motion mode changes cameras immediately.

The nine forms are a display gallery, not where these organisms live. Widths,
plaque thickness and each cell size are enlarged independently; cell sizes also
change between presets without biological meaning. Section patterns and the wavy
gum-epithelium boundary are artistic cues, not measured fibres, ridges or
trabeculae. The crevice is a separate authored explanatory patch, not measured
segmentation of the molar. The 2D fallback disables this 3D scale and retains
Pocket. Known limits are listed in docs/crevice-anatomy-verification.md.

Generate reusable local GLB objects with `node tools/export-crevice-objects.mjs
<output-directory>`. This exports the same nine forms with source identities,
normalized authoring units and hashes. They can be opened in Blender or another
GLB viewer. Exported forms never read saved density fields.
See docs/crevice-anatomy-verification.md for the current checks and limits;
docs/organic-crevice-verification.md and docs/crevice-verification.md record the
earlier passes.

Choose Species → Saved model to enlarge the map. Click a cell-form card, or
choose a taxon below the scene, to see its field and morphology source. “All”
shows dominant identity per bin; several taxa can coexist in a bin. Each
individual field has its own contrast normalization.

Evidence shows selected human specimen observations. Those positions describe
organization across attached plaque and do not create universal bands down a
pocket. Healthy does not inherit advanced-disease positions. Actinomyces remains
a group; spirochete observations do not identify every cell as T. denticola.

Mouth and Face use a BodyParts3D adult male reference atlas. The jaw opening is
an authored display pose. “Selected FDI 36” opens an independent authored tooth
specimen, not a measured registration. Tooth/Pocket Nerves & blood provides
schematic local routes and a named sensory pathway. Missing atlas nerve meshes
are not invented.

The optional biofilm story explains assembly and tissue-response concepts.
Its slider is stage position, not elapsed biological time. Bone loss is a
separate possible outcome, not inevitable progression from the earlier stages.

Drag to orbit, scroll to zoom, use arrow keys on the scene and Home to reset.
Keyboard focus is visible. Narrow layouts stack content and let the enlarged detail flow down the page. The 2D fallback is pocket.html?render=2d. The motion=reduce
query checks reduced-motion presentation without changing Mac settings.

Build with npm run build; check with npm run check. The existing organic pocket
build also requires Blender; runtime viewing uses bundled assets. Optional
source links open the internet. No private data are uploaded by this viewer.

Do not edit source hashes, original saved fields, solver kernels or raw atlas
vertices to improve the appearance. Work on presentation and source-separated
records. See docs/anatomy-map-verification.md for evidence and remaining limits.
Source control includes the viewer upgrade and its rebuild rules; no push is performed.
Generated atlas and face buffers remain local under the repository size policy.
See docs/reference-assets.md for their pinned sources and rebuild instructions.
