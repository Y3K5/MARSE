# Protein relations in a spatial lattice

This optional MARSE viewer reads an existing, retained protein-relation build. It displays source-coordinate protein backbones, literal evidence relationships, separate reference-identity links and a highlighted review queue. Independent network and protein cameras keep the surrounding network available while inspecting a fold or supported region. It never runs a scientific method.

## Inputs → what happens → outputs

Input is a local `marse-protein-relations-build-v1` folder produced by the upstream Workbench build interface. The importer checks every listed file, retains the verified JSON bytes without a second read, and refuses changed hashes, unsafe paths, mismatched source/model/trace identities and differing existing runtime data. Supplied JavaScript is never installed or executed.

The browser rechecks its imported producer manifest and caches the verified JSON bytes, so it does not reread unchecked replacements. The viewer checks its source/residue/packet pins, literal endpoints, reference identities and explicit residue-name/chain/number/insertion joins before rendering. Missing folds remain markers; unsupported partner sites stay unresolved. Output is an interactive, schematic presentation and a local camera/filter save file. Upstream Workbench owns scientific evidence qualification, study revisions and recorded-run custody. A matching checksum is not biological validation or authenticated authorship.

## Run locally

From the MARSE repository root, first install/check the existing studio runtime libraries as described in the parent guide. No new dependency or asset download is added here.

```sh
python3 tools/organism-studio/spatial-lattice/setup_relations.py --from <retained-build-folder>
python3 tools/organism-studio/serve.py
```

Open `http://127.0.0.1:8988/spatial-lattice/index.html`. The parent studio has a **Protein relations** link. Add `?review=queue` to start with the imported review queue. Existing identical runtime imports are idempotent; differing imports are refused. Preserve the existing runtime and use a fresh destination when testing another input. The `--output` option is a local CLI test/import destination, never a browser-supplied filesystem path.

Runtime evidence goes only in ignored `runtime/`. Neither the build payloads, private links, local catalogue, models, data nor study receipts are part of the reusable code. Do not force-add runtime files. This folder contains code, tests and documentation only.

A fresh clone needs the parent studio's exact, separately held runtime libraries and a local retained build; there is no approved public asset/data bundle. With no valid import the viewer fails closed and shows an import message.

## See and switch the current protein

Now viewing keeps the current protein or source-pair name, reference, coordinate status and position in the visible list above the canvas. Choose Protein group and Protein, or use Previous/Next. Next starts at the first visible record when the group overview is open; navigation stops at the ends. Group overview returns to the field. Protein browser (Browse on phones) opens search, the full list and advanced filters. Coordinate-held records stay labelled and cannot open a substituted fold.

Inspect residues opens the existing single-protein lens. Switching to another retained fold updates that lens; switching to a coordinate hold closes the previous fold. Published source pairs have their own selector and Previous/Next controls, preserving literal endpoints and source material limits. A source endpoint is never labelled as a retained-group member through a reference connector.

Details and evidence (Evidence on phones) opens a readable drawer below the pinned navigation/network section. Browser and evidence drawers keep the current name and three-microbe overview visible; close them with their button or Escape. Narrow windows collapse the group/protein selectors under Switch protein or group while keeping current identity and step controls visible. Search retains the existing scope: review searches stay within the queue; forms searches can include coordinate-held catalogue records in the selected reference.

Saved views retain an explicitly visible selected protein with the existing camera/source pins; changed filters, unknown identities and conflicting frames are refused. Old saves without a current-protein field remain valid. Viewing/navigation never adopts a scientific target, changes evidence or launches analysis.

## Spacious exploration

The review queue opens as a separate protein field with larger retained folds and clearance between display envelopes. Volume, layers and orbits keep true XYZ variation; Protein spacing changes the gap while Fold size changes the displayed form. These are display transforms, never measured molecular distances. A gap marker never receives a substitute fold. Main names appear on hover or selection to leave the folds visible. Click a fold to focus it; Fit protein field restores the whole current field. Full names and exact references stay in the list/inspector.

The independent source-network overview stays above the field and close-up and remains pinned while the workspace scrolls, showing every organism and the literal source atlas with the candidate layer off. Drag/zoom that overview separately or use Fit overview; its glyphs summarize the graph. Explore source relations adds the original source endpoints to the main field for the existing evidence/region inspector. Choosing a retained-form shortcut returns to the larger protein field. No proximity or organism shortcut asserts a protein correspondence or interaction.

Organism shortcuts use imported metadata and the already retained forms. Folds, coordinate holds, source strains, ranks and original reference qualifications keep their own status. All references combines the existing form entries; it does not select new targets or create a reviewed cohort. High detail and residue beads are available for small fields; larger groups start balanced. `?view=review`, `?view=forms:<reference-id>` and `?view=all` can reopen known views; the existing `?review=queue` link still opens the queue.

Fold sizes span1–4.5 and spacing0.7–3. Saved views retain both the main and optional overview camera; old saves without an overview camera remain valid. A manual overview camera is preserved across layout resizing; Fit overview restores the full atlas. Historical immutable producer views keep their original renderer and layout; camera replay across different renderer versions is not independently qualified.

## All microbes together

Choose **All microbes + evidence** or open `?view=together&layout=layers`. Combined fields pack by actual display extents, and a thicker schematic backbone stroke keeps folds readable at wide zoom without altering the centerline coordinates. Every existing retained form/reference hold appears with the literal published-source endpoints in one spacious field. The full review queue keeps its highlights, including its missing-coordinate markers even when a forms-only filter would normally omit them. Bright source-pair links and endpoint rings identify named-protein direct assay/biochemical evidence across organisms; the evidence buttons focus each original pair and open its source inspector. Other cross-microbe records stay faint with their own context class. Candidate records stay off. The pinned independent overview keeps the full source atlas visible even during close-ups.

Source folds are separate entries from the retained results. Dashed blue reference connectors remain identity/navigation only. A source strain, construct or interaction is never transferred by species, label, shared fold or display proximity; source strains/material appear in the inspector. Existing tested constructs/reactive bonds keep their own gates. A source pair does not become a solved contact site, native docking pose or a confirmed link to a review protein. No evidence, input coordinates or historical build is rewritten by this shortcut.

## Large pairs and supported locations

The supported source-pair buttons open both literal folds as a large side-by-side display, with the whole three-microbe source network pinned above and a separate site lens beside them. Pair size changes their display magnification; it never docks or aligns them. Every retained residue remains source-bound. Source strains/material still apply, and unavailable folds remain markers. All microbial folds returns to the complete combined field.

The lens starts at the first qualified region on one endpoint. Choose Endpoint and Supported region, then Zoom supported location or zoom with the wheel/buttons. Gold highlights the supported source residues against faded retained context. An assayed construct is not a minimal contact site; a reactive bond is not a solved protein-pair interface. Unsupported partner sites stay unresolved with site zoom disabled. Optional layers remain separately controlled and off by default. The site, pair and whole-network cameras move independently. Narrow windows stack the large pair and site panels; Supported-site close-up jumps to the lens while the network overview stays pinned.

Add `&pair=<literal-record-id>` to a combined-view link to open a known supported pair directly. Unknown, candidate, context-only or non-protein pairs are refused. Save view retains the pair choice/magnification and existing main/overview cameras. Reopen restores those and opens the pair lens at its first qualified region; the independent site camera/selection is not part of the camera save contract. Old views remain valid. No original study, packet or producer build changes.

## Select highlighted amino acids

The site lens labels small qualified sites directly. Select a gold bead or its label to see the amino-acid name, distinct paper/construct number, canonical-reference number and structure-author chain/number/insertion code. Selection reveals the compact numbering/source summary. The selected residue card shows the evidence category, original source locator/method/material, source limitations and exact retained Workbench revision. Evidence and numbering jumps to that card without changing the pair or overview camera. Focus this residue zooms the independent site camera.

Large constructs use a few collision-managed labels and retain every highlighted residue in the searchable, paged Highlighted amino acids list. Search a name or a paper/canonical/author number, then select a row. The selected label remains available. Lines connect labels to the actual qualified C-alpha positions; labels and list entries never add coordinates. Context residues outside the current annotation receive source identity only and no site-evidence assignment. Changing endpoint/region clears the previous selection.

A residue inside a tested construct has construct-membership evidence, not an individually measured binding role. Reactive bonds, secondary cleavages, conditional reference projections and predictions retain their own labels and claim limits. Unsupported counterparts have no site labels or residue evidence. Source links open only when explicitly selected; labels do not fetch external sources or upload data.

## Explore

Choose clusters, layers or orbits. Search and filter retained results; queue searches remain within the queue. Gold halos show review membership, amber a label-only source navigation hold and teal a mapped source-supported pair. Source strain/material and scientific limits remain separate from that reference match. A similar label never creates a correspondence or interaction edge. Candidate records and optional reference/prediction/secondary-cut layers are off by default.

Select a protein to inspect its retained result and original holds. **Inspect fold beside the network** opens an independent camera. A source-coordinate bead retains its amino-acid name, chain, author residue and insertion code. No side chains, missing residues or molecular surfaces are reconstructed. Site overlays accept only the existing supported evidence categories; solved protein-pair interfaces and guessed counterpart residues are refused.

Each supported region links to its producer's exact local Workbench revision. The imported manifest describes the retained producer build; the current reusable viewer is separate presentation code and does not claim to be that original renderer snapshot. Return links are restricted to loopback HTTP and keep their imported revision identity. External source links require HTTP(S), are explicitly selected by the user and are not fetched automatically. The viewer has no upload endpoint or external API calls.

## Check with invented inputs

The fixture generator creates synthetic software data, never a biological example or new prediction. Keep its output outside the tracked source, or in ignored runtime storage.

```sh
python3 tools/organism-studio/spatial-lattice/fixture_relations.py <new-test-build-folder>
python3 tools/organism-studio/spatial-lattice/setup_relations.py --from <new-test-build-folder> --output <new-test-runtime-folder>
MARSE_RELATIONS_FIXTURE=<new-test-runtime-folder> node --test tools/organism-studio/spatial-lattice/relations.test.mjs tools/organism-studio/spatial-lattice/exploration.test.mjs tools/organism-studio/spatial-lattice/pair.test.mjs tools/organism-studio/spatial-lattice/site-residue.test.mjs tools/organism-studio/spatial-lattice/navigation.test.mjs tools/organism-studio/spatial-lattice/neighborhood.test.mjs
PYTHONPATH=tools/organism-studio/spatial-lattice python -m pytest tools/organism-studio/spatial-lattice/test_setup_relations.py
```

The Python suite covers checksums, paths/symlinks, source/model/trace drift, JSON-only import, replacement-read refusal, idempotency and preservation. Node tests cover metadata-driven groups, literal edges, review/hold classification, all three layouts, residue joins, site/partner refusals and saved filters. Repeat browser checks after visual changes. Numerical-engine tests and biological qualification are separate.

## Ownership and limits

`app.js` owns navigation and drawer state; `neighborhood-model.js` builds the selected protein’s evidence slice and display layout using existing classification/reference/site gates; `navigation-model.js` supplies visible-list identities, source-pair choices, bounded step navigation and saved-selection scope checks; `model.js` and `network-model.js` own input/display contracts and schematic layouts; `review-model.js` and `subtractive-model.js` own review navigation/filtering. `residue-model.js`, `fold-geometry.js`, `interaction-model.js`, `interaction-geometry.js` and `interaction-lens.js` reuse retained coordinate/region presentation. `runtime-loader.js` checks/caches the imported producer payloads in the browser. `setup_relations.py` owns local transport checks, not publication or scientific reanalysis.

The layout does not establish native placement, abundance, docking geometry, binding strength, folding dynamics, target suitability or clinical effect. A source label is not an identity match; no displayed connection means unresolved coverage, not biological absence. Existing numerical engines and the separate private research viewer are unchanged. Saved files pin the input snapshot; cross-renderer camera replay is not independently qualified.

The whole-source overview retains unassigned source/context markers in a neutral color; these do not create an additional organism label or a molecular identity.

## Protein neighborhood

Select a protein, then choose Protein neighborhood (Neighbors on phones). The main field contains the selected protein, any existing exact-reference source correspondences and their literal supported named-protein partners. The whole source atlas stays in its independent pinned overview. Reference-linked retained and source folds are separate; their source strain/material qualifications remain. No label-only, homolog, cell/class, candidate or phenotype evidence becomes a supported partner. An unmatched protein keeps its fold or coordinate hold with an unresolved explanation, plus explicit navigation to separate source neighborhoods.

Gold overlays consume the current validated site reader. Choose a qualified region or its bead to open the independent site lens, residue identity and original evidence/revision. Tested constructs, reactive bonds and optional secondary/conditional/predicted layers retain their own labels; a partner without a qualified region says site unresolved. Overlapping constructs are separate annotations, not independent contact residues. No displayed edge specifies a docked pose. Other recorded contexts remain in a separate expandable evidence section.

Previous/Next and the protein picker move through the current neighborhood without changing its center. Center neighborhood here or Explore on a partner makes that protein the new center. Return to protein field restores the field camera when available. Saved neighborhoods pin the existing source snapshot, center and optional selected member; unknown centers, selections outside the neighborhood and conflicting pair/frame state are refused. Old saves remain valid.

Connections and qualified regions use the existing side browser; choose Connections & sites (Connections on phones) when it is hidden. This leaves the fold viewport clear. Selecting a region closes the drawer and opens the independent site lens. The neighborhood browser and evidence drawers keep the pinned overview visible.

In neighborhood site inspection, relationship/endpoint/region controls start collapsed under Change relationship, endpoint or region, leaving the chosen fold and evidence more visible. Expand that row to change the qualified selection.

Return-to-field replays a captured camera only from the molecular display field. Entering a neighborhood from source-coordinate space returns to a fitted field, preserving the distinction between Angstrom coordinates and display units.

## Optional atoms, bonds and molecular context

Atoms & bonds opens a separate local atom inspector beside the network. At a qualified residue, Atoms & bonds here keeps its original paper/canonical/author identity. Choose a source model in this view, chain, residue and atom; click an atom or its bonded-neighbor list to follow topology. Local environment within6Å shows geometric surroundings, not additional assigned contacts. Fit whole model returns to the retained coordinate context. Source charges, occupancy, stored B/confidence fields and experimental-versus-predicted model status remain explicit; no missing hydrogens, solvent, protonation or native partner assembly is reconstructed.

The optional setup_atomic.py importer requires the existing Biopython structural parser. A private CLI configuration supplies pinned retained PDB view/mmCIF and component-template files plus the current runtime and a fresh output directory. It reads verified bytes once, checks every atom against its raw source and every C-alpha/name/chain/insertion anchor against this exact fold, and writes a separate immutable ignored runtime/atomic payload. It never changes original runtime data or native studies. No source JavaScript or new renderer is installed; the inspector uses the existing local Three runtime. Use a fresh successor output rather than overwriting any import.

Browser atomic loading checks the source snapshot, import/payload/PDB hashes and exact atom/backbone correspondence. No exact full-atom import means a visible hold; related proteins or reference connectors cannot supply a substitute. Component bond orders are retained chemical-component targets; polymer links are derived only from consecutive declared polymer positions and source geometry. Source disulfide, covalent, metal coordination and hydrogen connection records remain distinct. These are not atomic protein-pair interface evidence, docking or a chemical-reaction calculation. Atomic camera/selection is an independent temporary inspection and does not change the existing saved-view contract.


Atom ownership: `setup_atomic.py` qualifies the private import, `atomic-model.js` checks source bytes/joins and supplies topology queries, and `atomic-lens.js` renders the independently navigable close-up. `app.js` opens it from the selected protein or supported residue and clears its state when leaving. The atom inspector does not create a StudyStore revision or native MARSE scientific build.

Run `node --test tools/organism-studio/spatial-lattice/atomic.test.mjs` for invented atom/name/element/chain/insertion/coordinate, bond, stale-pin, safe-path and no-substitution checks. With the existing structural-parser environment, run `PYTHONPATH=tools/organism-studio/spatial-lattice python -m pytest tools/organism-studio/spatial-lattice/test_setup_atomic.py`. Import with `python tools/organism-studio/spatial-lattice/setup_atomic.py --config <private-config.json>`. Configuration schema is `marse-atomic-import-config-v1`: `runtime`, fresh `output`, `template_path`/`template_sha256`, and `models` containing `protein_id`, `view_path`/`view_sha256`, `raw_path`/`raw_sha256` and public source `context` (`id`, `kind`, `source_url`, method/title/scope/admission). Keep this configuration and all resulting files local and ignored. No configuration belongs in a saved browser view.

Format definitions: [chemical-component bond orders](https://mmcif.wwpdb.org/dictionaries/mmcif_ma.dic/Items/_chem_comp_bond.value_order.html), [source connection types](https://mmcif.wwpdb.org/dictionaries/mmcif_ma.dic/Items/_struct_conn.conn_type_id.html), and [polymer sequence identity](https://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_atom_site.label_seq_id.html). Author numbering remains distinct from label sequence positions and any qualified paper/canonical numbering.

Checked locally: 4 atom-reader groups, 7 importer/topology tests, 1 renderer geometry/disposal group, 33 existing viewer groups and 13 importer controls. Actual retained imports and desktop/phone display checks have private receipts with their research owner. These checks qualify the software and retained-source correspondence, not the biological validity of a predicted fold or native atom-level interface. Full-atom coverage is optional per exact fold; the rest of the catalogue remains backbone-only or held. No new scientific computation or dependency installation was performed.


Commit review hardened incomplete numeric atom records and experimental structure identity, preserved the neighborhood camera when changing only the atomic model, and verified instance-buffer disposal through repeated selections. Fresh imports describe retained water from their actual atom records and refuse existing/dangling-symlink output locations. Existing imported evidence is preserved. `node --test tools/organism-studio/spatial-lattice/atomic-geometry.test.mjs` checks geometry and buffer disposal using the already-installed local Three assets. The core package dependencies are unchanged; the optional importer uses a separately installed Biopython parser.
