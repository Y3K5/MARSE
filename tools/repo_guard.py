#!/usr/bin/env python3
"""Structural and supply-chain guard for the MARSE repository.

Branch protection and access control live in GitHub's settings, where only the
repository owner can set them (see GOVERNANCE.md). This script enforces the
things settings cannot: that the repository's *contents* stay safe and
structured, on every commit and in CI.

It checks two families of problem.

**Workflow security.** A GitHub Actions workflow is code that runs with access
to the repository and its secrets, so a careless or hostile change to one is
the most valuable target in a repository like this. The guard requires every
action to be pinned to a full commit SHA, refuses the ``pull_request_target``
trigger, refuses write permissions that were not deliberately declared, refuses
untrusted interpolation into shell commands, and refuses any publishing step
that is not gated behind a GitHub environment.

**Structure.** Required files must exist, the package layout must match the
documented architecture, and unexpected files must not accumulate at the
repository root.

Usage, from the repository root::

    python tools/repo_guard.py

It exits 0 when everything passes and 1 with a list of problems otherwise.
Only the standard library is used, so it runs anywhere Python does.
"""

from __future__ import annotations

import re
import sys
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

# --- workflow rules ---------------------------------------------------------

PINNED_ACTION = re.compile(r"^[\w.-]+/[\w.-]+(?:/[\w.-]+)*@[0-9a-f]{40}$")
"""owner/repo@<40 hex> - a tag or branch can be moved, a commit SHA cannot."""

USES_LINE = re.compile(r"^\s*-?\s*uses:\s*(?P<ref>\S+)")
PERMISSION_LINE = re.compile(r"^\s*(?P<name>[a-z-]+):\s*(?P<value>read|write|none)\s*$")

FORBIDDEN_TRIGGERS = {
    "pull_request_target": (
        "runs with the base repository's secrets while checking out a fork's code; "
        "use 'pull_request' instead"
    ),
}

PUBLISHING_MARKERS = (
    "pypa/gh-action-pypi-publish",
    "softprops/action-gh-release",
    "ncipollo/release-action",
    "twine upload",
    "gh release create",
    "poetry publish",
    "flit publish",
    "npm publish",
)
"""Steps that publish outside the repository. Each needs an approval gate."""

# Untrusted text interpolated straight into a shell command is a script
# injection: a branch or issue title can carry shell metacharacters.
UNTRUSTED_INTERPOLATION = re.compile(
    r"\$\{\{\s*(github\.event\b|github\.head_ref\b|inputs\.)[^}]*\}\}"
)

ALLOWED_WRITE_PERMISSIONS: dict[str, set[str]] = {
    # workflow file -> permissions it may legitimately hold as write
}
"""Deliberate exceptions. Empty: no workflow currently needs write access."""

# --- structure rules --------------------------------------------------------

REQUIRED_FILES = (
    "README.md",
    "LICENSE",
    "PRIVACY.md",
    "SECURITY.md",
    "GOVERNANCE.md",
    "CONTRIBUTING.md",
    "CHANGELOG.md",
    "CITATION.cff",
    "pyproject.toml",
    ".gitignore",
    ".gitattributes",
    ".pre-commit-config.yaml",
    ".github/CODEOWNERS",
    ".github/workflows/ci.yml",
    ".github/workflows/privacy.yml",
    "docs/specification.md",
    "docs/architecture.md",
    "docs/theory.md",
    "docs/parameters.md",
    "docs/validation.md",
    "docs/roadmap.md",
    "tools/privacy_guard.py",
    "tools/repo_guard.py",
)

REQUIRED_PACKAGE_DIRS = (
    "core",
    "schemas",
    "spatial",
    "microbes",
    "biofilm",
    "adaptation",
    "interventions",
    "validation",
    "ecosystem",
    "analysis",
    "evidence",
    "experimental",
)
"""The architecture from docs/architecture.md. Removing one is a design change."""

ALLOWED_ROOT_ENTRIES = frozenset(
    {
        ".git",
        ".github",
        ".gitattributes",
        ".gitignore",
        ".pre-commit-config.yaml",
        ".privacy-allowlist",
        "CHANGELOG.md",
        "CITATION.cff",
        "CONTRIBUTING.md",
        "GOVERNANCE.md",
        "LICENSE",
        "PRIVACY.md",
        "README.md",
        "SECURITY.md",
        "docs",
        "examples",
        "pyproject.toml",
        "src",
        "tests",
        "tools",
    }
)
"""Tracked entries permitted at the repository root, to stop it sprawling."""


@dataclass(frozen=True)
class Problem:
    where: str
    detail: str

    def __str__(self) -> str:
        return f"{self.where}: {self.detail}"


def _strip_comment(line: str) -> str:
    """Remove a trailing YAML comment, ignoring '#' inside quotes."""
    out, quote = [], ""
    for char in line:
        if quote:
            out.append(char)
            if char == quote:
                quote = ""
        elif char in "\"'":
            quote = char
            out.append(char)
        elif char == "#":
            break
        else:
            out.append(char)
    return "".join(out).rstrip()


def _indent(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def _workflow_files(root: Path) -> list[Path]:
    directory = root / ".github" / "workflows"
    if not directory.is_dir():
        return []
    return sorted(p for p in directory.iterdir() if p.suffix in {".yml", ".yaml"})


def _permission_blocks(lines: list[str]) -> Iterator[tuple[int, str, str]]:
    """Yield (line number, permission name, value) for every permissions: block."""
    for index, raw in enumerate(lines):
        line = _strip_comment(raw)
        if not re.match(r"^\s*permissions:\s*$", line):
            continue
        base = _indent(line)
        for offset in range(index + 1, len(lines)):
            entry = _strip_comment(lines[offset])
            if not entry.strip():
                continue
            if _indent(entry) <= base:
                break
            match = PERMISSION_LINE.match(entry)
            if match:
                yield offset + 1, match["name"], match["value"]


def _shell_blocks(lines: list[str]) -> Iterator[tuple[int, str]]:
    """Yield (line number, body) for each ``run:`` step, body only.

    A block scalar (``run: |``) runs until the indentation returns to the
    ``run:`` key's level, so neighbouring keys such as a ``with:`` input on the
    next step are not part of it.
    """
    pattern = re.compile(r"^(?P<lead>\s*(?:-\s+)?)(?:run|script):(?P<inline>.*)$")
    for index, raw in enumerate(lines):
        match = pattern.match(_strip_comment(raw))
        if not match:
            continue
        inline = match["inline"].strip()
        if inline and inline not in {"|", ">", "|-", ">-", "|+", ">+"}:
            yield index + 1, inline
            continue
        base = len(match["lead"])
        body: list[str] = []
        for offset in range(index + 1, len(lines)):
            following = lines[offset]
            if following.strip() and _indent(following) <= base:
                break
            body.append(following)
        yield index + 1, "\n".join(body)


def check_workflow(path: Path, relative: str) -> list[Problem]:
    problems: list[Problem] = []
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()

    # 1. Every action pinned to an immutable commit.
    for number, raw in enumerate(lines, start=1):
        match = USES_LINE.match(_strip_comment(raw))
        if not match:
            continue
        ref = match["ref"].strip("\"'")
        if ref.startswith("./") or ref.startswith("docker://"):
            continue
        if not PINNED_ACTION.match(ref):
            problems.append(
                Problem(
                    f"{relative}:{number}",
                    f"action '{ref}' is not pinned to a 40-character commit SHA; "
                    "a tag can be moved under you",
                )
            )

    # 2. Triggers that expose secrets to untrusted code.
    for number, raw in enumerate(lines, start=1):
        line = _strip_comment(raw)
        for trigger, reason in FORBIDDEN_TRIGGERS.items():
            if re.match(rf"^\s*{trigger}\s*:", line):
                problems.append(Problem(f"{relative}:{number}", f"'{trigger}' {reason}"))

    # 3. Write permissions must be declared as deliberate exceptions.
    allowed = ALLOWED_WRITE_PERMISSIONS.get(relative, set())
    has_permissions = False
    for number, name, value in _permission_blocks(lines):
        has_permissions = True
        if value == "write" and name not in allowed:
            problems.append(
                Problem(
                    f"{relative}:{number}",
                    f"grants '{name}: write'; add it to ALLOWED_WRITE_PERMISSIONS in "
                    "tools/repo_guard.py if that is intended",
                )
            )
    if not has_permissions:
        problems.append(
            Problem(relative, "no 'permissions:' block; declare least privilege explicitly")
        )

    # 4. Untrusted text interpolated into a shell command. Only the run block's own
    #    body counts: an expression in a neighbouring 'with:' input is not shell.
    for number, body in _shell_blocks(lines):
        found = UNTRUSTED_INTERPOLATION.search(body)
        if found:
            problems.append(
                Problem(
                    f"{relative}:{number}",
                    f"interpolates {found.group()} into a shell command; "
                    "pass it through an 'env:' variable and quote it instead",
                )
            )

    # 5. Publishing must be gated behind an environment that can require approval.
    lowered = text.lower()
    for marker in PUBLISHING_MARKERS:
        if marker in lowered:
            if not re.search(r"^\s*environment:", text, flags=re.M):
                problems.append(
                    Problem(
                        relative,
                        f"publishes ('{marker}') without an 'environment:'; a protected "
                        "environment is what lets GitHub require approval before release",
                    )
                )
            break
    return problems


def check_workflows(root: Path) -> list[Problem]:
    problems: list[Problem] = []
    workflows = _workflow_files(root)
    if not workflows:
        return [Problem(".github/workflows", "no workflows found; CI protects this repository")]
    for path in workflows:
        problems += check_workflow(path, path.relative_to(root).as_posix())
    return problems


def check_structure(root: Path, tracked: Iterable[str]) -> list[Problem]:
    problems: list[Problem] = []

    for required in REQUIRED_FILES:
        if not (root / required).is_file():
            problems.append(Problem(required, "required file is missing"))

    package = root / "src" / "marse"
    if not package.is_dir():
        problems.append(Problem("src/marse", "package directory is missing"))
    else:
        for name in REQUIRED_PACKAGE_DIRS:
            subpackage = package / name
            if not (subpackage / "__init__.py").is_file():
                problems.append(
                    Problem(
                        f"src/marse/{name}",
                        "documented architecture subpackage is missing its __init__.py",
                    )
                )
        for module in sorted(package.rglob("*.py")):
            if module.name == "py.typed":
                continue
            text = module.read_text(encoding="utf-8").lstrip()
            if not text.startswith(('"""', "'''")):
                problems.append(
                    Problem(
                        module.relative_to(root).as_posix(),
                        "module has no docstring; every module states its responsibility",
                    )
                )

    roots = {entry.split("/", 1)[0] for entry in tracked}
    for unexpected in sorted(roots - ALLOWED_ROOT_ENTRIES):
        problems.append(
            Problem(
                unexpected,
                "unexpected entry at the repository root; keep new work inside an existing "
                "directory, or add it to ALLOWED_ROOT_ENTRIES in tools/repo_guard.py",
            )
        )
    return problems


def tracked_files(root: Path) -> list[str]:
    """Files git tracks, or every file on disk when git is unavailable."""
    import subprocess

    try:
        completed = subprocess.run(
            ["git", "ls-files", "-z"], cwd=root, capture_output=True, check=True
        )
    except (OSError, subprocess.CalledProcessError):
        return [
            p.relative_to(root).as_posix()
            for p in root.rglob("*")
            if p.is_file() and ".git" not in p.parts
        ]
    return [entry for entry in completed.stdout.decode().split("\0") if entry]


def main(argv: list[str] | None = None) -> int:
    root = Path(argv[0]).resolve() if argv else Path.cwd()
    problems = check_workflows(root) + check_structure(root, tracked_files(root))
    if problems:
        for problem in problems:
            print(problem)
        print(f"\nrepo guard: {len(problems)} problem(s). GOVERNANCE.md explains the rules.")
        return 1
    workflows = len(_workflow_files(root))
    print(f"repo guard: {workflows} workflow(s) and the repository structure are sound")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
