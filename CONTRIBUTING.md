# Contributing to MARSE

Thank you for helping build MARSE. The project is in Phase 0 (specification),
so design discussion in issues is as valuable as code right now.

## Ground rules

- **Discuss first.** Open an issue, using a template, before a large change and
  before any new scientific model.
- **Respect the scope.** MARSE v1.0 is deliberately narrow
  ([specification](docs/specification.md)). Ideas beyond it are welcome as
  issues for later.
- **Protect privacy.** Read [PRIVACY.md](PRIVACY.md) and complete its one-time
  setup before your first commit. Commits must use a GitHub noreply email
  address; the hooks and CI enforce this.

## Development setup

Requires Python 3.12 or newer.

```bash
git clone https://github.com/Y3K5/MARSE.git    # or your fork
cd MARSE
python -m venv .venv
source .venv/bin/activate                       # Windows: .venv\Scripts\activate
python -m pip install -e ".[dev]"
pre-commit install
```

## Workflow

1. Create a branch from `main`: `git switch -c short-description`.
2. Make focused commits with clear messages.
3. Run `pre-commit run --all-files` and `python -m pytest`.
4. Update the documentation and the *Unreleased* section of `CHANGELOG.md`.
5. Open a pull request and complete its checklist.

CI runs the pre-commit hooks and the tests on Linux (Python 3.12 to 3.14) and
Windows, and the privacy workflow runs on every push.

## Code standards

- Formatting and linting use ruff, through pre-commit.
- Public functions have type hints and NumPy-style docstrings.
- The core stays deterministic: randomness only from seeded generators owned by
  the core, no global random state, no network access, and nothing
  machine-specific written into outputs.
- A new runtime dependency is justified in its pull request; the core's
  dependencies stay minimal.

## Tests

Run `python -m pytest`. Markers mirror the [validation plan](docs/validation.md):

| Marker | Purpose |
|---|---|
| *(none)* | Unit tests: each component implements its stated rule |
| `numerical` | Agreement with analytical or high-precision references |
| `invariance` | Results unchanged under problem-preserving transformations |
| `regression` | Comparison with stored golden trajectories |
| `benchmark` | Documented biological benchmark cases (slow) |

For example, `python -m pytest -m "not benchmark"` skips the slow cases.

## Scientific contributions

A new model, rule or parameter set states:

- its biological scale and its spatial and temporal resolution;
- its equations or rules, with units;
- its assumptions and known limitations;
- a published source and a confidence note for every parameter;
- at least one validation case. Benchmark definitions come before new biology.

## Review and release

Every change reaches the default branch through a pull request approved by the
maintainer; nobody merges their own work. Releases are a deliberate, gated act.
See [GOVERNANCE.md](GOVERNANCE.md).

## Maintenance

- GitHub Actions are pinned to commit SHAs and updated by Dependabot.
- Pre-commit hooks are pinned the same way; update them with
  `pre-commit autoupdate --freeze`. When gitleaks changes version, also update
  `GITLEAKS_VERSION` and `GITLEAKS_SHA256` in
  [privacy.yml](.github/workflows/privacy.yml) from that release's checksums file.
- Before a release, follow [PRIVACY.md](PRIVACY.md#publishing).

## License

Contributions are licensed under the [Apache License 2.0](LICENSE), as its
section 5 describes.
