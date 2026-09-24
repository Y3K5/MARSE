# Governance

Who decides what happens to MARSE, who can change it, and how a release is
made. The short version: **one maintainer approves every change, and nothing
is published without a deliberate, gated act.**

## Roles

| Role | Who | Can |
|---|---|---|
| Maintainer | [@Y3K5](https://github.com/Y3K5) | Merge to the default branch, tag releases, change repository settings, grant access |
| Contributor | Anyone | Open issues, open pull requests from a fork |

There is currently one maintainer. Adding another is a deliberate decision
recorded by a change to this file, not an ad-hoc settings change.

## Changes

1. Work happens on a branch or a fork, never directly on the default branch.
2. Every change reaches the default branch through a pull request.
3. Every pull request needs the maintainer's approval
   ([`.github/CODEOWNERS`](.github/CODEOWNERS) requests it automatically).
4. CI and the privacy workflow must pass first.

Contributors cannot merge their own work, and nobody can bypass review by
pushing straight to the default branch, once the settings below are in place.

## Releasing

**No workflow in this repository publishes anything, and none may be added
without an approval gate.** [`tools/repo_guard.py`](tools/repo_guard.py) fails
CI if a publishing step appears without a GitHub `environment:`, which is what
allows GitHub to require a human approval before the step runs.

A release is therefore an explicit act by the maintainer:

1. Confirm the privacy and publication checklist in
   [PRIVACY.md](PRIVACY.md#publishing) — this is the point at which names and
   contact details become permanently public.
2. Update `CHANGELOG.md` and the version in `src/marse/__init__.py`.
3. Tag the release, and build artifacts in CI from a clean checkout.
4. Archive the release (for example on Zenodo) only when the content is final.
   Archived records are permanent.

Nothing is published to a package index yet. When that changes, it goes
through a protected environment with the maintainer as a required reviewer,
and trusted publishing rather than a long-lived API token.

## What the repository enforces itself

Run on every commit through pre-commit, and again in CI:

| Guard | Enforces |
|---|---|
| [`tools/repo_guard.py`](tools/repo_guard.py) | Actions pinned to commit SHAs; no `pull_request_target`; no undeclared write permissions; no untrusted interpolation into shell commands; no ungated publishing step; required files present; documented package layout intact; no sprawl at the repository root |
| [`tools/privacy_guard.py`](tools/privacy_guard.py) | No personal data, credentials, blocked file types or oversized files; commit identities use a noreply address |
| gitleaks | No secrets, including across the full history in CI |

These protect the repository's *contents*. They cannot stop someone who
already has write access, which is what the settings below are for.

## Settings only the owner can apply

GitHub enforces these, not the repository. They are listed here so the
intended posture is written down and auditable, and they are the part of
[PRIVACY.md's checklist](PRIVACY.md#github-settings-checklist) that concerns
access rather than disclosure.
[docs/repository-setup.md](docs/repository-setup.md) walks through applying
them, including which are unavailable on which plan.

**Access**

- [ ] Two-factor authentication on the owner account. An account takeover
      defeats everything else on this page.
- [ ] No collaborators beyond the maintainer; no organisation-wide write access.
- [ ] Deploy keys: none, or read-only.

**Default branch ruleset** (Settings → Rules → Rulesets, targeting the default branch)

> Rulesets and classic branch protection are **not available for private
> repositories on GitHub Free**; they need Pro, Team or Enterprise Cloud. They
> are free for *public* repositories on every plan. While this repository is
> private on a free plan, the checks below cannot be enforced by GitHub, and
> the repository's own guards plus the maintainer's discipline are what stand
> in for them. See [docs/repository-setup.md](docs/repository-setup.md).

- [ ] Require a pull request before merging.
- [ ] Require review from Code Owners.
- [ ] Require the `CI` and `Privacy` status checks to pass.
- [ ] Block force pushes.
- [ ] Restrict deletion.
- [ ] Do not allow bypassing the above, including for the owner. Bypass
      permissions are how protections quietly stop applying.

**Actions** (Settings → Actions → General)

- [ ] Workflow permissions: read-only by default.
- [ ] Do not allow GitHub Actions to create or approve pull requests.
- [ ] Require approval for all outside collaborators' workflow runs.

**Publishing**

- [ ] Restrict tag creation with a **tag ruleset** targeting `v*`. The old tag
      protection rules were sunset in 2024 and migrated into rulesets, so the
      same plan limitation applies.
- [ ] Create the protected environment before the first publishing workflow
      exists, not after.

**Security** (Settings → Code security)

- [ ] Dependabot alerts: on. Available on every plan.
- [ ] Private vulnerability reporting: **public repositories only.** Enable it
      as part of going public; until then there is no external reporter, since
      nobody else can see the repository.
- [ ] Secret scanning with push protection: free on public repositories. On
      private repositories it needs the paid GitHub Secret Protection product,
      so gitleaks in CI is what covers this repository meanwhile — it scans the
      full history on every push and needs no licence.

**Visibility**

- [ ] The repository is private while it is being built. Going public is a
      one-way step: forks and clones survive a later switch back. Complete the
      checklist in [PRIVACY.md](PRIVACY.md#github-settings-checklist) first.

## Reporting a problem

Security and privacy problems go through
[private reporting](https://github.com/Y3K5/MARSE/security/advisories/new), not
a public issue. See [SECURITY.md](SECURITY.md).

## Changing this document

Amending governance is itself a pull request, reviewed like any other change.
