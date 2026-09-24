# Privacy

MARSE is developed in the open, so this repository rests on one assumption:
**anything that reaches GitHub is public forever.** Deleting a file later does
not remove it from git history, forks, clones, search engines, software
archives such as Software Heritage, or releases archived on Zenodo. Protection
therefore happens *before* anything leaves your machine, and CI checks again
as a second line of defence.

## What never enters the repository

| Never commit | Why | Instead |
|---|---|---|
| Personal email address, phone number, postal address | Published in files and in commit metadata | Your GitHub noreply address; contact through GitHub |
| Your name, until you choose to publish it | You decide when to be identified | Your GitHub handle; add names at release time ([Publishing](#publishing)) |
| Passwords, tokens, API keys, SSH or GPG keys, `.env` files | Anyone could use them | Environment variables and a password manager |
| Private, unpublished or human-derived data | Research ethics, participant privacy, priority of your results | Keep it outside git and document how to obtain it ([Data](#data)) |
| Personal notes, plans and drafts | Internal thinking | The git-ignored `private/` folder |
| PDF and Office documents | Hidden author names, usernames and edit history | Markdown; figures as PNG or SVG |
| Photos with EXIF or XMP metadata | GPS location, device serial number, owner name | Strip the metadata first |
| Absolute paths such as `/Users/<name>/…` or `C:\Users\<name>\…` | They reveal your username | Paths relative to the repository or experiment file |
| Notebook outputs, log files, `.DS_Store` | They capture paths, usernames, data and file names | Clear outputs; keep logs local |

## How this is enforced

| Layer | When | What it checks |
|---|---|---|
| 1. Git hooks on your machine | Before each commit | Secrets (gitleaks); personal emails, home-folder paths, blocked file types, oversized files, notebook outputs and image metadata (privacy guard); your private denylist; the identity git will record |
| | When you write a commit message | Emails, paths and denylist entries in the message |
| | Before each push | The author and committer identity of every outgoing commit |
| 2. Your GitHub account | When you push | GitHub rejects pushes that expose your private email address |
| 3. CI: [privacy.yml](.github/workflows/privacy.yml) | After every push, on every branch | The same checks across all files and every commit, plus a gitleaks scan of the full history |
| 4. Repository settings | Always | Secret scanning with push protection, private vulnerability reporting ([checklist](#github-settings-checklist)) |

Layers 1 and 2 prevent leaks. Layer 3 is a backstop: on a public repository,
anything it catches has already been published, so treat a failure there as an
incident ([If something leaks anyway](#if-something-leaks-anyway)).

The privacy guard is [tools/privacy_guard.py](tools/privacy_guard.py). Its
output masks email addresses and never repeats a denylist entry, so it is safe
in public CI logs.

## One-time setup

1. **GitHub account.** Under *Settings → Emails*, turn on **Keep my email
   addresses private** and **Block command line pushes that expose my email**.
   The same page shows your noreply address; for this account it is
   `38446806+Y3K5@users.noreply.github.com`.

2. **Git identity** in your clone of this repository:

   ```bash
   git config user.name "Y3K5"
   git config user.email "38446806+Y3K5@users.noreply.github.com"
   ```

   Add `--global` to use this identity for every repository.

3. **Install the hooks:**

   ```bash
   python -m pip install -e ".[dev]"
   pre-commit install
   ```

   This installs the pre-commit, commit-msg and pre-push hooks. The first run
   downloads the tools; gitleaks is built with Go, which pre-commit fetches if
   you do not have it.

4. **Private denylist** (recommended). Create a file *outside* the repository
   that lists exact strings that must never be committed: your full name and
   its variants, personal email addresses, phone numbers, your home address,
   institutional login names.

   - macOS and Linux: `~/.config/marse/private-denylist.txt`
   - Windows: `%USERPROFILE%\.config\marse\private-denylist.txt`

   Put one entry per line; `#` starts a comment and matching ignores case. To
   keep the file elsewhere, set `MARSE_PRIVATE_DENYLIST_FILE` to its path. The
   guard reports only *where* an entry matched, never the entry itself.

5. **Denylist in CI** (optional). Add the same lines as a repository secret
   named `PRIVATE_DENYLIST` (*Settings → Secrets and variables → Actions*).
   Secrets are not passed to pull requests from forks, so those runs skip this
   one check.

## The private folder

`private/` at the repository root is ignored by git. Keep plans, drafts, notes
and documents such as a development plan PDF there. If a file inside it is
ever force-added, the guard refuses it.

## Data

- Keep datasets outside the repository, or in the git-ignored `data/` folder,
  and document how to obtain them.
- Commit only data that is small and either synthetic or already published
  under a license that allows redistribution, for example golden regression
  trajectories produced by MARSE itself.
- Never commit human-derived data or metadata about people: clinical isolates
  with patient information, microbiome samples containing host reads,
  participant records.
- Parameter values cite published sources, not unpublished lab notebooks.

## Privacy by design in MARSE itself

Run manifests are meant to be shared with papers, issues and archives, so MARSE
must never write identifying information into them:

- no usernames, hostnames, IP addresses or environment variables;
- no absolute paths: inputs are recorded relative to the experiment file,
  together with a SHA-256 checksum;
- platform details limited to the operating-system family and the versions of
  Python and key libraries;
- timestamps in UTC, so the local time zone is not revealed.

MARSE makes no network requests and collects no telemetry. These rules are part
of the [architecture](docs/architecture.md#privacy-rules-for-manifests-and-outputs)
and will be enforced by tests.

## When a check blocks you

| Message | What to do |
|---|---|
| `... is not a GitHub noreply address` | Set your identity (setup step 2). To fix unpushed commits: `git rebase -r origin/main --exec "git commit --amend --no-edit --reset-author"` |
| `personal email address` or `matches ... private denylist` | Remove it; refer to people by their GitHub handle |
| `absolute home-directory path` | Use a path relative to the repository or experiment file |
| `never commit this kind of file` | Move it to `private/` or outside the repository. If it truly belongs, strip its metadata and list it in [`.privacy-allowlist`](.privacy-allowlist) |
| `image has EXIF metadata` or `XMP metadata` | `exiftool -all= -overwrite_original <image>`, or re-export the image (plotting libraries and screenshots usually add none) |
| `notebook has saved outputs` | `jupyter nbconvert --clear-output --inplace <notebook>` |
| `file is ... MiB` | Keep the data outside git; commit a script or download instructions instead |
| A gitleaks finding | Remove the secret **and revoke it** with its provider: assume it is compromised |

Do not bypass the hooks with `--no-verify`.

## If something leaks anyway

1. **Revoke secrets first.** Rotate the key or token with its provider straight
   away. Rewriting history does not un-leak it.
2. **Remove it from history** with [git filter-repo](https://github.com/newren/git-filter-repo)
   (for example `--invert-paths --path <file>` or `--replace-text`), then
   force-push every affected branch and tag.
3. **Ask GitHub Support** to purge cached views and pull-request references,
   which survive a force-push.
4. **Check forks and releases.** Zenodo records are permanent; contact Zenodo
   support if a release archive is affected.
5. **Assume it was copied.** Removal reduces exposure but cannot guarantee
   deletion, so also change what you can (a new email address, new keys).

## Publishing

Publication is a deliberate, one-way disclosure. Before the first tagged
release or paper submission, decide what you want to make public:

- Names, ORCID iDs and affiliations go into `CITATION.cff` and later
  `paper/paper.md`. An ORCID iD credits you without exposing an email address.
- If a contact email is required, use an institutional or dedicated project
  address, not a personal inbox.
- Build release artifacts in CI from a clean checkout. The source distribution
  uses an explicit file list, so untracked files in a working copy are never
  shipped.
- Run `pre-commit run --all-files` and review every identity-related change
  before tagging.

## GitHub settings checklist

Account:

- [ ] Emails: keep addresses private and block pushes that expose them.
- [ ] Public profile: review your name, location, company, bio and linked accounts.
- [ ] Two-factor authentication: on. An account takeover would expose everything.
- [ ] Commit signing, if you use it: an SSH signing key, or a GPG key whose only
      identity is your noreply address. Public GPG keys list their email addresses.

Repository settings:

- [ ] Security: turn on Dependabot alerts. Private vulnerability reporting and
      free secret scanning arrive when the repository goes public; gitleaks in
      CI covers secret scanning until then.
- [ ] Actions: workflow permissions set to read-only, and approval required for
      workflows from outside contributors' forks.
- [ ] Rules: protect `main` by requiring the CI and Privacy checks and blocking
      force-pushes. Needs a paid plan while the repository is private, and is
      free once it is public — see
      [docs/repository-setup.md](docs/repository-setup.md).
- [ ] Features: turn off the wiki unless you need it; wiki edits bypass all of
      these checks.
- [ ] Secrets: `PRIVATE_DENYLIST` (optional, setup step 5).
- [ ] Make the repository public only after everything above is done and the
      Privacy workflow is green.
