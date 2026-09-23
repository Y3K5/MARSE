#!/usr/bin/env python3
"""Privacy guard for the MARSE repository.

Keeps personal information and private research material out of git history.
It runs locally through pre-commit, before anything leaves your machine, and
again in CI as a backstop. Secrets such as API keys and tokens are detected by
gitleaks, which runs alongside this script.

Usage, from the repository root::

    python tools/privacy_guard.py files [--all] [PATH ...]
    python tools/privacy_guard.py identity
    python tools/privacy_guard.py message COMMIT_MSG_FILE
    python tools/privacy_guard.py commits [--outgoing] [REV ...]

Output never repeats private denylist entries and masks email addresses, so
it is safe to show in public CI logs. PRIVACY.md describes the policy this
script enforces. Only the standard library is used.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import re
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

MAX_FILE_BYTES = 1024 * 1024
"""Larger files are refused: datasets and bulk outputs belong outside git."""

MAX_SCAN_BYTES = 8 * 1024 * 1024
"""Upper bound on how much of any single file is read."""

DENYLIST_ENV = "MARSE_PRIVATE_DENYLIST"
DENYLIST_FILE_ENV = "MARSE_PRIVATE_DENYLIST_FILE"
ALLOWLIST_FILE = ".privacy-allowlist"
PRIVATE_DIR = "private/"

# Commit identities must use one of these addresses; none of them reveals an inbox.
NOREPLY_EMAIL = re.compile(
    r"""^(?:
        [^@\s]+@users\.noreply\.github\.com   # GitHub's per-account noreply address, bots
      | noreply@github\.com                   # commits made in the GitHub web interface
      | noreply@anthropic\.com                # commits made by the Claude coding assistant
    )$""",
    re.IGNORECASE | re.VERBOSE,
)

# Addresses that may appear in file contents in addition to the noreply ones.
PLACEHOLDER_EMAIL = re.compile(
    r"""^(?:
        git@github\.com                                    # SSH clone URLs
      | [^@\s]+@(?:[a-z0-9-]+\.)*example\.(?:com|net|org)  # RFC 2606 documentation domains
      | [^@\s]+@(?:[a-z0-9-]+\.)*(?:example|invalid|localhost|test)
    )$""",
    re.IGNORECASE | re.VERBOSE,
)

EMAIL = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)*\.[A-Za-z]{2,}")

# "logo@2x.png" looks like an address; skip matches whose "TLD" is a file suffix.
FILE_SUFFIX_TLDS = frozenset(
    {"csv", "gif", "jpeg", "jpg", "json", "md", "pdf", "png", "py", "svg", "toml", "txt", "webp"}
    | {"yaml", "yml"}
)

HOME_PATH = re.compile(
    r"(?<![\w.-])(?:/Users/|/home/|[A-Za-z]:\\{1,2}Users\\{1,2})([A-Za-z0-9][A-Za-z0-9._-]*)",
    re.IGNORECASE,
)
GENERIC_HOME_NAMES = frozenset(
    {"admin", "codespace", "default", "example", "jovyan", "me", "name", "public", "root"}
    | {"runner", "shared", "ubuntu", "user", "username", "vscode", "you", "yourname"}
    | {"your-name", "your_name"}
)

# Why each kind of file is refused, then its file-name patterns (matched case-insensitively).
BLOCKED_FILE_TYPES: tuple[tuple[str, str], ...] = (
    (
        "credential or key file",
        ".env .env.* .envrc .netrc _netrc .pypirc .git-credentials credentials.json id_rsa* "
        "id_dsa* id_ecdsa* id_ed25519* *.pem *.key *.p12 *.pfx *.jks *.keystore *.kdbx *.gpg "
        "*.asc *.ovpn *private-denylist*",
    ),
    (
        "PDF or office document (these carry hidden author metadata)",
        "*.pdf *.doc *.docx *.docm *.xls *.xlsx *.xlsm *.ppt *.pptx *.odt *.ods *.odp *.rtf "
        "*.pages *.numbers",
    ),
    (
        "raw research data (datasets stay outside git)",
        "*.fastq *.fq *.fastq.gz *.fq.gz *.sam *.bam *.bai *.cram *.crai *.vcf *.vcf.gz *.bcf "
        "*.sra *.fast5 *.pod5 *.czi *.nd2 *.lif *.lsm",
    ),
    (
        "archive or database (its contents cannot be reviewed)",
        "*.zip *.tar *.tgz *.gz *.bz2 *.xz *.7z *.rar *.db *.sqlite *.sqlite3",
    ),
    (
        "OS or log file (reveals usernames, paths or file names)",
        ".DS_Store Thumbs.db desktop.ini *.log",
    ),
)
TEMPLATE_FILE_NAMES = frozenset({".env.example", ".env.sample", ".env.template"})

EXIF_MARKERS: dict[str, tuple[bytes, ...]] = {
    ".jpg": (b"Exif\x00\x00",),
    ".jpeg": (b"Exif\x00\x00",),
    ".png": (b"eXIf", b"Raw profile type exif"),
    ".webp": (b"EXIF",),
    ".heic": (b"Exif",),
    ".heif": (b"Exif",),
    ".avif": (b"Exif",),
}
XMP_START = b"adobe:ns:meta/"
XMP_END = b"xmpmeta>"
XMP_SENSITIVE_FIELDS = (b"GPS", b"dc:creator", b"Artist", b"Author", b"Owner", b"SerialNumber")

LOG_FORMAT = "%H%x1f%an%x1f%ae%x1f%cn%x1f%ce%x1e"
IDENT = re.compile(r"^(?P<name>.*) <(?P<email>[^>]*)> \d+ [+-]\d{4}$")


class GitError(RuntimeError):
    """A git command failed."""


@dataclass(frozen=True)
class Finding:
    where: str
    problem: str

    def __str__(self) -> str:
        return f"{self.where}: {self.problem}"


def git(*args: str, cwd: Path | None = None) -> str:
    completed = subprocess.run(["git", *args], cwd=cwd, capture_output=True, check=False)
    if completed.returncode != 0:
        raise GitError(completed.stderr.decode("utf-8", "replace").strip())
    return completed.stdout.decode("utf-8", "replace")


def repo_root() -> Path:
    return Path(git("rev-parse", "--show-toplevel").strip()).resolve()


def mask_email(address: str) -> str:
    """Hide most of ``address`` so that reports never re-publish it."""
    local, _, domain = address.partition("@")
    host, dot, tld = domain.rpartition(".")
    return f"{local[:1]}***@{host[:1]}***{dot}{tld}"


def redact(text: str, denylist: Sequence[str]) -> str:
    for term in denylist:
        text = re.sub(re.escape(term), "***", text, flags=re.IGNORECASE)
    return text


def matches_denylist(text: str, denylist: Sequence[str]) -> bool:
    folded = text.casefold()
    return any(term in folded for term in denylist)


def denylist_path(environ: Mapping[str, str]) -> Path | None:
    if environ.get(DENYLIST_FILE_ENV):
        return Path(environ[DENYLIST_FILE_ENV]).expanduser()
    config_home = environ.get("XDG_CONFIG_HOME")
    if not config_home:
        try:
            config_home = str(Path.home() / ".config")
        except RuntimeError:  # no home directory to look in
            return None
    return Path(config_home) / "marse" / "private-denylist.txt"


def load_denylist(environ: Mapping[str, str]) -> list[str]:
    """Personal terms (names, emails, phone numbers...) that must never be committed.

    They come from a file outside the repository and, in CI, from a secret.
    """
    sources = [environ.get(DENYLIST_ENV, "")]
    path = denylist_path(environ)
    if path is not None and path.is_file():
        sources.append(path.read_text(encoding="utf-8-sig"))
    terms = {line.strip().casefold() for source in sources for line in source.splitlines()}
    return sorted(term for term in terms if term and not term.startswith("#"))


def load_allowlist(root: Path) -> list[str]:
    path = root / ALLOWLIST_FILE
    if not path.is_file():
        return []
    lines = (line.strip() for line in path.read_text(encoding="utf-8").splitlines())
    return [line for line in lines if line and not line.startswith("#")]


def is_publishable_email(address: str) -> bool:
    return bool(NOREPLY_EMAIL.match(address) or PLACEHOLDER_EMAIL.match(address))


def blocked_reason(rel_path: str) -> str | None:
    """Why the file at ``rel_path`` ('/'-separated, from the root) may never be committed."""
    lowered = rel_path.lower()
    if lowered.startswith(PRIVATE_DIR):
        return "file inside the git-ignored private/ folder"
    name = lowered.rsplit("/", 1)[-1]
    if name in TEMPLATE_FILE_NAMES:
        return None
    for reason, patterns in BLOCKED_FILE_TYPES:
        if any(fnmatch.fnmatchcase(name, pattern.lower()) for pattern in patterns.split()):
            return reason
    return None


def scan_text(
    text: str, where: str, denylist: Sequence[str], *, generic: bool = True
) -> list[Finding]:
    """Look for personal email addresses, home-directory paths and denylisted terms."""
    check_denylist = bool(denylist) and matches_denylist(text, denylist)
    findings: list[Finding] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        location = f"{where}:{lineno}"
        if generic and "@" in line:
            for match in EMAIL.finditer(line):
                address = match.group()
                if address.rpartition(".")[2].lower() in FILE_SUFFIX_TLDS:
                    continue
                if not is_publishable_email(address):
                    problem = f"personal email address {mask_email(address)}"
                    findings.append(Finding(location, problem))
        if generic:
            for match in HOME_PATH.finditer(line):
                if match.group(1).lower() not in GENERIC_HOME_NAMES:
                    problem = (
                        "absolute home-directory path (reveals a username); use a relative path"
                    )
                    findings.append(Finding(location, problem))
        if check_denylist and matches_denylist(line, denylist):
            findings.append(Finding(location, "matches an entry in your private denylist"))
    return findings


def notebook_problem(text: str) -> str | None:
    try:
        cells = json.loads(text).get("cells", [])
    except (ValueError, AttributeError):
        return None
    if any(isinstance(cell, dict) and cell.get("outputs") for cell in cells):
        return (
            "notebook has saved outputs (they often capture paths, usernames or data); clear "
            "them with: jupyter nbconvert --clear-output --inplace <notebook>"
        )
    return None


def image_metadata_problems(data: bytes, suffix: str) -> list[str]:
    problems = []
    if any(marker in data for marker in EXIF_MARKERS.get(suffix, ())):
        problems.append(
            "image has EXIF metadata (can include GPS location, device serial number, owner); "
            "strip it before committing"
        )
    start = data.find(XMP_START)
    if suffix in EXIF_MARKERS and start != -1:
        end = data.find(XMP_END, start)
        packet = data[start : end if end != -1 else len(data)]
        if any(field in packet for field in XMP_SENSITIVE_FIELDS):
            problems.append(
                "image has XMP metadata with location or authorship fields; strip it "
                "before committing"
            )
    if suffix == ".svg" and b"<dc:creator" in data:
        problems.append("SVG metadata names its creator; remove the <metadata> block")
    return problems


def check_file(
    root: Path, rel: str, allowlist: Sequence[str], denylist: Sequence[str]
) -> list[Finding]:
    """Check one file, given by its '/'-separated path relative to ``root``."""
    findings: list[Finding] = []
    if denylist and matches_denylist(rel, denylist):
        findings.append(Finding(rel, "file name matches an entry in your private denylist"))
    exempt = any(fnmatch.fnmatchcase(rel, pattern) for pattern in allowlist)
    reason = None if exempt else blocked_reason(rel)
    if reason:
        return [*findings, Finding(rel, f"never commit this kind of file: {reason}")]

    path = root / rel
    if path.is_symlink():
        target = os.readlink(path)
        return findings + scan_text(target, f"{rel} (symlink target)", denylist, generic=not exempt)
    if not path.is_file():
        return findings

    size = path.stat().st_size
    if size > MAX_FILE_BYTES and not exempt:
        findings.append(
            Finding(
                rel,
                f"file is {size / 1_048_576:.1f} MiB (limit {MAX_FILE_BYTES / 1_048_576:g} MiB); "
                "keep datasets and bulk outputs outside git",
            )
        )
    with path.open("rb") as handle:
        data = handle.read(MAX_SCAN_BYTES)
    suffix = path.suffix.lower()
    if not exempt:
        findings += [Finding(rel, problem) for problem in image_metadata_problems(data, suffix)]

    text = data.decode("utf-8", errors="replace")
    if b"\0" in data[:8192]:  # binary: only the denylist applies to its bytes
        if denylist and matches_denylist(text, denylist):
            findings.append(
                Finding(rel, "binary file contains an entry from your private denylist")
            )
        return findings
    findings += scan_text(text, rel, denylist, generic=not exempt)
    if suffix == ".ipynb" and not exempt:
        problem = notebook_problem(text)
        if problem:
            findings.append(Finding(rel, problem))
    return findings


def relative_to_root(path: str, root: Path) -> str:
    candidate = Path(path)
    # Resolve only the parent: resolving a symlink itself would follow it.
    absolute = candidate.parent.resolve() / candidate.name
    try:
        return absolute.relative_to(root).as_posix()
    except ValueError:
        return candidate.as_posix()


def run_files(paths: Sequence[str], all_files: bool, denylist: Sequence[str]) -> list[Finding]:
    root = repo_root()
    if all_files:
        rels = [rel for rel in git("ls-files", "-z", cwd=root).split("\0") if rel]
        status = "active" if denylist else "not configured"
        print(f"privacy guard: checking {len(rels)} tracked files (private denylist: {status})")
    else:
        rels = [relative_to_root(path, root) for path in paths]
    allowlist = load_allowlist(root)
    findings: list[Finding] = []
    for rel in rels:
        findings += check_file(root, rel, allowlist, denylist)
    return findings


def identity_findings(
    where: str, role: str, name: str, email: str, denylist: Sequence[str]
) -> list[Finding]:
    findings = []
    if not NOREPLY_EMAIL.match(email):
        findings.append(
            Finding(
                where,
                f"{role} email {mask_email(email)} is not a GitHub noreply address, and git "
                "publishes it with every commit (see PRIVACY.md)",
            )
        )
    if denylist and matches_denylist(f"{name}\n{email}", denylist):
        findings.append(Finding(where, f"{role} name or email matches your private denylist"))
    return findings


def run_identity(denylist: Sequence[str]) -> list[Finding]:
    """Check the identity git will record on the next commit."""
    findings: list[Finding] = []
    for role in ("author", "committer"):
        try:
            ident = git("var", f"GIT_{role.upper()}_IDENT").strip()
        except GitError:
            problem = f"no {role} identity configured; set user.name and user.email (PRIVACY.md)"
            findings.append(Finding("git config", problem))
            continue
        match = IDENT.match(ident)
        name, email = (match["name"], match["email"]) if match else (ident, "")
        findings += identity_findings("git config", role, name, email, denylist)
    return findings


def outgoing_revs(environ: Mapping[str, str]) -> list[str]:
    from_ref = environ.get("PRE_COMMIT_FROM_REF")
    to_ref = environ.get("PRE_COMMIT_TO_REF")
    if from_ref and to_ref:
        return [f"{from_ref}..{to_ref}"]
    return ["--branches", "--not", "--remotes"]  # everything not yet on a remote


def run_commits(revs: Sequence[str], denylist: Sequence[str]) -> list[Finding]:
    if not revs:
        try:
            git("rev-parse", "--verify", "--quiet", "HEAD")
        except GitError:  # no commits yet
            return []
        revs = ["HEAD"]
    log = git("log", f"--format={LOG_FORMAT}", *revs, "--")
    findings: list[Finding] = []
    for record in log.split("\x1e"):
        fields = record.strip("\n").split("\x1f")
        if len(fields) != 5:
            continue
        sha, author, author_email, committer, committer_email = fields
        where = f"commit {sha[:12]}"
        findings += identity_findings(where, "author", author, author_email, denylist)
        findings += identity_findings(where, "committer", committer, committer_email, denylist)
    return findings


def run_message(path: str, denylist: Sequence[str]) -> list[Finding]:
    kept = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("# ") and ">8" in line:  # scissors line; the diff below it is not sent
            break
        kept.append("" if line.startswith("#") else line)
    return scan_text("\n".join(kept), "commit message", denylist)


def report(findings: Sequence[Finding], denylist: Sequence[str]) -> int:
    if not findings:
        return 0
    for finding in findings:
        print(redact(str(finding), denylist))
    print(f"\nprivacy guard: {len(findings)} problem(s). PRIVACY.md explains how to fix each kind.")
    return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="privacy_guard.py",
        description="Keep personal information and private research material out of git.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    files = commands.add_parser("files", help="check file names, sizes and contents")
    files.add_argument("paths", nargs="*", help="files to check")
    files.add_argument("--all", action="store_true", help="check every file tracked by git")
    commands.add_parser("identity", help="check the identity git will record on the next commit")
    message = commands.add_parser("message", help="check a commit message file")
    message.add_argument("file")
    commits = commands.add_parser("commits", help="check commit author and committer identities")
    commits.add_argument("revs", nargs="*", help="revisions for git log (default: HEAD)")
    commits.add_argument(
        "--outgoing", action="store_true", help="check the commits about to be pushed"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    args = build_parser().parse_args(argv)
    denylist = load_denylist(os.environ)
    try:
        if args.command == "files":
            findings = run_files(args.paths, args.all, denylist)
        elif args.command == "identity":
            findings = run_identity(denylist)
        elif args.command == "message":
            findings = run_message(args.file, denylist)
        else:
            revs = outgoing_revs(os.environ) if args.outgoing else args.revs
            findings = run_commits(revs, denylist)
    except GitError as error:
        print(f"privacy guard: git failed: {redact(str(error), denylist)}", file=sys.stderr)
        return 2
    return report(findings, denylist)


if __name__ == "__main__":
    sys.exit(main())
