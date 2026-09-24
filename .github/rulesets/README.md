# Rulesets

Branch and tag protection for this repository, kept as JSON so that a change to
who can write to the default branch is reviewable like any other change.

GitHub enforces these; the repository cannot. See
[GOVERNANCE.md](../../GOVERNANCE.md) for the whole picture and
[docs/repository-setup.md](../../docs/repository-setup.md) for the other
settings.

## Importing

**Settings → Rules → Rulesets → New ruleset ▾ → Import a ruleset**, then pick
the file and click **Create**. Do this once per file:

| File | Protects |
|---|---|
| `protect-default-branch.json` | The default branch: pull request required, checks must pass, no force push, no deletion |
| `protect-release-tags.json` | `v*` tags: cannot be moved or deleted once cut |

Rulesets are free on public repositories. They are unavailable for private
repositories on GitHub Free.

## Why zero required approvals

`required_approving_review_count` is **0**, and `require_code_owner_review` is
**false**. That looks weaker than it should be, and it is deliberate:

**GitHub does not let anyone approve their own pull request.** With a single
maintainer, requiring one approval — or requiring code-owner review, since the
code owner is also the only author — means no pull request can ever be merged.
Combined with an empty bypass list, that locks the maintainer out of the
default branch entirely.

What the ruleset still buys with zero approvals is the part that matters most
here: **every change reaches the default branch through a pull request, with
CI and the privacy checks green**. Nothing lands unchecked, and nothing can be
force-pushed or deleted.

### When a second maintainer joins

Change two values in `protect-default-branch.json` and re-import:

```json
"required_approving_review_count": 1,
"require_code_owner_review": true
```

At that point [CODEOWNERS](../CODEOWNERS) becomes binding and self-merging
stops being possible, which is the behaviour the governance document
describes for a multi-maintainer project.

## The bypass list

**Ruleset exports deliberately omit the bypass list**, so importing these files
gives you no bypass at all. That is the intended state: a bypass entry for the
maintainer is the most common reason a ruleset silently stops protecting
anything.

If you ever need one, add it in the UI after importing, and treat it as
temporary.

## Other settings that matter

`strict_required_status_checks_policy` is **false**, so a pull request need not
be rebased onto the latest default branch before merging. With one maintainer
the risk of a semantic conflict is small, and `true` forces a CI re-run after
every unrelated merge. Set it to `true` once several people are merging in
parallel.

The seven required checks are the job names CI publishes. If a workflow job is
renamed, the requirement silently stops matching — rename it here too, and
confirm in a pull request that the check still appears as required.

## Verifying

After importing, push a trivial commit directly to the default branch. GitHub
should refuse it. If it succeeds, check that the ruleset's enforcement is
**Active** rather than *Evaluate*, and that its bypass list is empty.
