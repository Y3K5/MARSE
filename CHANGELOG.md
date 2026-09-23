# Changelog

All notable changes to MARSE are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and versions follow
[Semantic Versioning](https://semver.org/). Any change that alters simulation
results for the same manifest is always called out.

## [Unreleased]

### Added

- Package layout following the planned architecture, and a `marse --version`
  command.
- Scientific specification, architecture, validation plan and roadmap in `docs/`.
- Privacy protection: local git hooks (gitleaks, plus a privacy guard for files,
  commit messages, git identity and outgoing commits) and a CI workflow that
  checks every file and every commit and scans the full history for secrets.
- CI running the pre-commit hooks and the tests on Linux (Python 3.12 to 3.14)
  and Windows.
- Apache-2.0 license and citation metadata.
