# Code-only publication candidate

Runtime assets are excluded through the studio .gitignore: vendor/, structures/, integration/ and community.json. ASSET_MANIFEST.json pins exact paths, sizes and SHA-256 hashes. setup_assets.py imports an existing local asset set without network access and refuses missing or modified files. serve.py checks assets before starting. Original source bytes, embedded public depositor emails and license notices remain intact locally.

The candidate publishes code, tests, source references, documentation and the asset manifest only. No generated engine output, raw coordinate file, third-party runtime bundle, illustration, historical private report or machine-specific script is selected. The private studio remains preserved. No privacy allowlist was relaxed.

A fresh clone needs the exact local assets before running. There is no automatic public download or approved asset redistribution bundle. Full third-party bundle notice verification and owner release decisions remain prerequisites for any later asset distribution; they do not claim approval through this code-only change.

Eight catalogue and six continuous checks passed with local assets before splitting. Final split verification: 76 pinned local assets, eight catalogue tests, six continuous groups and the setup import/refusal test pass. The repository privacy gate passes for the exact code-only file list. gitleaks is not installed and its separate secret scan remains unperformed; these are software checks, not biological qualification. Native placement, dynamics and calibration remain unqualified. The numerical engine is unchanged.

This branch starts from the existing local mac-visualizer HEAD; current remote-main compatibility is not established. No commit or push is authorized by this preparation. Review the exact candidate inventory before committing; the user performs every push.

## PR 69 validation update

The published head passed GitHub secrets and personal-data checks. GitHub and a local merge-tree check report no conflict with the fetched main revision. CI lint findings in the Python helpers and two final newlines were corrected locally; the new commit still requires a push and fresh CI before merge.
