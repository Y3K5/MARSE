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

## Explore

Choose clusters, layers or orbits. Search and filter retained results; queue searches remain within the queue. Gold halos show review membership, amber a label-only source navigation hold and teal a mapped source-supported pair. Source strain/material and scientific limits remain separate from that reference match. A similar label never creates a correspondence or interaction edge. Candidate records and optional reference/prediction/secondary-cut layers are off by default.

Select a protein to inspect its retained result and original holds. **Inspect fold beside the network** opens an independent camera. A source-coordinate bead retains its amino-acid name, chain, author residue and insertion code. No side chains, missing residues or molecular surfaces are reconstructed. Site overlays accept only the existing supported evidence categories; solved protein-pair interfaces and guessed counterpart residues are refused.

Each supported region links to its producer's exact local Workbench revision. The imported manifest describes the retained producer build; the current reusable viewer is separate presentation code and does not claim to be that original renderer snapshot. Return links are restricted to loopback HTTP and keep their imported revision identity. External source links require HTTP(S), are explicitly selected by the user and are not fetched automatically. The viewer has no upload endpoint or external API calls.

## Check with invented inputs

The fixture generator creates synthetic software data, never a biological example or new prediction. Keep its output outside the tracked source, or in ignored runtime storage.

```sh
python3 tools/organism-studio/spatial-lattice/fixture_relations.py <new-test-build-folder>
python3 tools/organism-studio/spatial-lattice/setup_relations.py --from <new-test-build-folder> --output <new-test-runtime-folder>
MARSE_RELATIONS_FIXTURE=<new-test-runtime-folder> node --test tools/organism-studio/spatial-lattice/relations.test.mjs
PYTHONPATH=tools/organism-studio/spatial-lattice python -m pytest tools/organism-studio/spatial-lattice/test_setup_relations.py
```

The Python suite covers checksums, paths/symlinks, source/model/trace drift, JSON-only import, replacement-read refusal, idempotency and preservation. Node tests cover metadata-driven groups, literal edges, review/hold classification, all three layouts, residue joins, site/partner refusals and saved filters. Repeat browser checks after visual changes. Numerical-engine tests and biological qualification are separate.

## Ownership and limits

`app.js` owns navigation; `model.js` and `network-model.js` own input/display contracts and schematic layouts; `review-model.js` and `subtractive-model.js` own review navigation/filtering. `residue-model.js`, `fold-geometry.js`, `interaction-model.js`, `interaction-geometry.js` and `interaction-lens.js` reuse retained coordinate/region presentation. `runtime-loader.js` checks/caches the imported producer payloads in the browser. `setup_relations.py` owns local transport checks, not publication or scientific reanalysis.

The layout does not establish native placement, abundance, docking geometry, binding strength, folding dynamics, target suitability or clinical effect. A source label is not an identity match; no displayed connection means unresolved coverage, not biological absence. Existing numerical engines and the separate private research viewer are unchanged. Saved files pin the input snapshot; cross-renderer camera replay is not independently qualified.
