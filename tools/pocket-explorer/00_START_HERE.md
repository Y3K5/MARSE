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
http://127.0.0.1:8975/pocket.html?view=crevice directly. **Opened section** exposes
the walls and attachment; **Look down** views the entrance; **Biofilm surface**
cuts away opposing tissue to inspect nine selectable 3D taxon/group examples.
Click a form or use the species selector to focus it and open its morphology
reference. The crevice now bends through its depth and tapers into the support.
Its gingival margin is a continuous fold. Exposed surfaces use softer highlights;
deliberate section faces are matte. Camera buttons travel gently for about half a
second; dragging, scrolling or arrow keys stop the travel. Reduced-motion mode
changes cameras immediately. Biofilm surface opens the lining and core while
keeping margin, ligament and bone context. The tissue patterns are artistic
surface cues, not measured fibers or trabeculae. Other forms remain present. These are a display library, not a
measured natural arrangement. Widths and each cell size are independently
exaggerated; sizes also change to fit the preset, without biological meaning.
The crevice is a separate authored explanatory patch, not measured segmentation
of the molar. The 2D fallback disables this 3D scale and retains Pocket.

Generate reusable local GLB objects with `node tools/export-crevice-objects.mjs
<output-directory>`. This exports the same nine forms with source identities,
normalized authoring units and hashes. They can be opened in Blender or another
GLB viewer. Exported forms never read saved density fields.
See docs/organic-crevice-verification.md for the current checks and limits;
docs/crevice-verification.md records the preceding scaffold pass.

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
