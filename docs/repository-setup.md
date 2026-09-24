# Repository setup runbook

How to apply the protections listed in [GOVERNANCE.md](../GOVERNANCE.md).
Written for the maintainer doing it by hand, once.

## The mental model

Repository protection splits in two, and only one half can live in the code:

| | Enforced by | Costs | Status |
|---|---|---|---|
| **What the contents allow** | `tools/repo_guard.py`, `tools/privacy_guard.py`, gitleaks, CI | Nothing | **Done.** Runs on every commit and in CI |
| **Who can do what** | GitHub settings | Plan-dependent | **Yours to apply** |

The first half is why this repository can refuse an unpinned action, an
ungated publishing step or a leaked address without anyone watching. It cannot
stop somebody who already has write access — that is what the second half is
for.

## Check your plan first

Most of the access controls are **not available for private repositories on
GitHub Free**. They need Pro, Team or Enterprise Cloud. The same features are
**free for public repositories on every plan**.

Check at **github.com → your avatar → Settings → Billing and plans**.

| Feature | Free + private | Free + public | Pro/Team + private |
|---|---|---|---|
| Two-factor authentication | yes | yes | yes |
| Actions permission settings | yes | yes | yes |
| Dependabot alerts | yes | yes | yes |
| Branch and tag rulesets | **no** | yes | yes |
| Secret scanning, push protection | **paid add-on** | yes | paid add-on |
| Private vulnerability reporting | **no** | yes | no |

Two consequences worth absorbing:

1. **On a free private repository, most of this is simply unavailable.** That
   is not a reason to skip the rest; it is a reason to know which guarantees
   you actually have. Right now the repository's own guards and your own
   discipline are standing in for the rulesets.
2. **Going public makes nearly all of it free**, and the project intends to be
   public anyway. Publishing is the cheapest way to get these protections, not
   a trade-off against them — provided the disclosure checklist in
   [PRIVACY.md](../PRIVACY.md#publishing) is done first.

---

## 1. Two-factor authentication

Do this one whatever your plan. An account takeover defeats every other
control on this page, because the attacker becomes you.

**Settings → Password and authentication → Two-factor authentication → Enable.**

Prefer an authenticator app or a passkey over SMS, which is vulnerable to SIM
swapping. **Save the recovery codes somewhere you can reach without GitHub.**
Losing them and your device means losing the account.

## 2. Actions permissions

Available on every plan, and worth setting before any workflow is added that
might want more.

**Repository → Settings → Actions → General:**

- *Workflow permissions* → **Read repository contents and packages
  permissions**. This is the default token scope; a workflow that genuinely
  needs more must ask for it in its own file, where `repo_guard.py` will see
  it and fail unless the exception is declared.
- Untick **Allow GitHub Actions to create and approve pull requests.** A
  workflow approving its own pull request would defeat required review.
- *Fork pull request workflows* → require approval for **all outside
  collaborators**, so a fork's workflow cannot run until you have read it.

## 3. Dependabot alerts

**Settings → Code security → Dependabot alerts → Enable.**

Available on every plan. This is the alerting half; the update half is already
configured in [`.github/dependabot.yml`](../.github/dependabot.yml), which
keeps the SHA-pinned actions current.

## 4. The default-branch ruleset

*Needs a paid plan while private, free once public.* This is the rule that
makes [`CODEOWNERS`](../.github/CODEOWNERS) binding rather than advisory —
without it, the file only **requests** review.

The configuration is committed as JSON, so import it rather than clicking
through it:

**Settings → Rules → Rulesets → New ruleset ▾ → Import a ruleset**, choose
[`.github/rulesets/protect-default-branch.json`](../.github/rulesets/protect-default-branch.json),
then **Create**.

That ruleset requires a pull request with passing checks, and blocks force
pushes and deletion of the default branch. It sets **zero required approvals
on purpose**: GitHub does not allow anyone to approve their own pull request,
so requiring one approval as a lone maintainer would lock you out of your own
default branch.
[`.github/rulesets/README.md`](../.github/rulesets/README.md) explains that in
full, and gives the two values to change when a second maintainer joins.

Then verify it actually bites: push a trivial commit straight to the default
branch and confirm GitHub refuses it. If it succeeds, check that enforcement
is **Active** rather than *Evaluate*, and that the bypass list is empty.

> **A one-person repository still benefits from this.** Not because you do not
> trust yourself, but because it removes the class of mistake where a tired
> `git push` to the wrong branch bypasses every check the project relies on.

## 5. Tag ruleset for releases

*Same plan requirement.* Tag protection rules were sunset in 2024 and migrated
into rulesets, so a tag ruleset is now the way to restrict who can cut a
release.

Import [`.github/rulesets/protect-release-tags.json`](../.github/rulesets/protect-release-tags.json)
the same way. It stops a `v*` tag being moved or deleted once cut, so a
published release cannot be quietly repointed at different code.

## 6. When the repository goes public

Do these in order, and only after the disclosure checklist in
[PRIVACY.md](../PRIVACY.md#publishing):

1. **Settings → General → Danger Zone → Change visibility → Public.**
   This is effectively one-way: forks, clones and caches survive a later
   switch back to private.
2. **Settings → Code security**, now free:
   - **Private vulnerability reporting → Enable.** This is what makes the
     reporting link in [SECURITY.md](../SECURITY.md) work.
   - **Secret scanning → Enable**, and **push protection → Enable**. This adds
     a check at push time; gitleaks in CI keeps scanning the full history.
3. Apply sections 4 and 5 above if you could not before, since they become
   free at this point.

## What to do while the rulesets are unavailable

If you are on a free plan and staying private for now, the honest position is
that GitHub is enforcing very little and the project is relying on:

- the repository's own guards, which run before every commit and in CI;
- there being exactly one person with write access;
- two-factor authentication on that person's account.

That is a reasonable posture for a solo project that is not yet public. It
stops being reasonable the moment a second person gets write access, or the
first external contribution arrives — at which point the ruleset is not
optional, and that is a good moment to go public and get it for free.

## Verifying it worked

| Check | Expected |
|---|---|
| Push directly to the default branch | Refused once the ruleset is active |
| Open a pull request | Review requested from `@Y3K5` automatically |
| Merge with a failing check | Blocked |
| `python tools/repo_guard.py` | Passes |
| `pre-commit run --all-files` | All hooks pass |
