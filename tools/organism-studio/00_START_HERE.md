# Experimental MARSE organism studio

Select an organism → inspect a connected surface patch → inspect source-linked protein coordinates in one Three.js scene. Drag to orbit and choose Whole cell to return. Source hashes, chains and construct caveats remain visible. Missing coordinates remain unavailable.

## Run locally

Runtime assets are deliberately Git-ignored. First run `python3 tools/organism-studio/setup_assets.py --from <existing-studio-folder>` using an existing local studio containing the exact asset versions. The command verifies every input before copying and makes no network requests. Run it without `--from` to check the installed assets. A fresh clone alone is not a turnkey distribution; no public asset bundle has been approved.

From the repository root: `python3 tools/organism-studio/serve.py`. Open http://127.0.0.1:8988/unified.html . Stop with Control-C. No installation or remote API is required. This loopback server serves only the studio directory.

## Check

From this directory run `node --test catalog.test.mjs` and `node check-continuous.mjs .`, plus `python3 test_setup_assets.py`. These verify data contracts, source availability and continuous geometry behavior. Browser inspection is still required after visual changes.

## Development

app.js owns selection; flow-dive.js owns the connected surface, source model and depth camera; source-protein.js verifies and parses local structures. unified.html/css/js provide navigation. Vendored Three.js provides rendering; 3Dmol supplies parsing in the continuous path. The numerical Python engine is separate. Community reads stored model fields or a user-selected local frame file; it does not run an engine simulation.

## Limits

Experimental visualization only. Cell geometry, membrane spacing, lipids, magnification and the linked inspection cradle are illustrative. Native insertion, sidedness, placement, abundance and binding are not established. Source backbone interpolation is not a secondary-structure assignment. Atomic display bonds are parser-derived aids. Saved model fields are exploratory, not measured human data or calibrated biology. No molecular dynamics is implemented.

The public candidate excludes historical local reports, machine-specific scripts, prior UI designs and generated illustrations. The private working original is preserved. See THIRD_PARTY.md and PUBLICATION_REVIEW.md before any release.
