# Changelog

All notable changes to MARSE are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and versions follow
[Semantic Versioning](https://semver.org/). Any change that alters simulation
results for the same manifest is always called out.

## [Unreleased]

### Added

- Growth kinetics: Monod, Haldane–Andrews substrate inhibition, Pirt uptake
  with maintenance, Luedeking–Piret product formation, and the Baranyi–Roberts
  growth curve with lag, evaluated in a numerically stable form.
- Secondary models: cardinal temperature (CTMI) and pH (CPM) models and the
  Ratkowsky square-root model, following the gamma concept.
- Oxygen solubility (Benson–Krause) and diffusivity (Han–Bartels) correlations,
  which refuse to extrapolate outside their stated validity ranges.
- Analytical reference solutions for validation cases V1–V5, each checked in
  the test suite against an independent numerical solution.
- `docs/theory.md` (every equation, its assumptions, numerical methods and
  limitations) and `docs/parameters.md` (values with units, sources and
  explicit confidence levels).
- Runnable examples: `batch_growth.py` and `oxygen_penetration.py`.
- Package layout following the planned architecture, and a `marse --version`
  command.
- Scientific specification, architecture, validation plan and roadmap in `docs/`.
- Privacy protection: local git hooks (gitleaks, plus a privacy guard for files,
  commit messages, git identity and outgoing commits) and a CI workflow that
  checks every file and every commit and scans the full history for secrets.
- CI running the pre-commit hooks and the tests on Linux (Python 3.12 to 3.14)
  and Windows.
- Apache-2.0 license and citation metadata.
