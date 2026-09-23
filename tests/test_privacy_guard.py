"""Tests for tools/privacy_guard.py, the repository's privacy gate.

Sample personal data is assembled at runtime (``addr("jane", "mail.zz")``,
``home("jane")``) so that this file never contains a literal address or home
path and passes the guard it tests. All names and numbers are fictional.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NOREPLY = "38446806+Y3K5@users.noreply.github.com"
SECRET_NAME = "Jane Q. Example"


def _load_guard():
    spec = importlib.util.spec_from_file_location(
        "privacy_guard", ROOT / "tools" / "privacy_guard.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


guard = _load_guard()


def addr(local: str, domain: str) -> str:
    return f"{local}@{domain}"


def home(*parts: str) -> str:
    return "/".join(["", "Users", *parts])


def git(cwd: Path, *args: str) -> str:
    completed = subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True, text=True)
    return completed.stdout.strip()


def commit(repo: Path, email: str, message: str = "change", name: str = "Y3K5") -> str:
    identity = ("-c", f"user.name={name}", "-c", f"user.email={email}")
    git(repo, *identity, "commit", "-q", "--allow-empty", "-m", message)
    return git(repo, "rev-parse", "HEAD")


@pytest.fixture(autouse=True)
def isolated_environment(tmp_path, monkeypatch):
    """Keep the developer's real denylist and git configuration out of the tests."""
    monkeypatch.setenv(guard.DENYLIST_FILE_ENV, str(tmp_path / "no-denylist.txt"))
    monkeypatch.delenv(guard.DENYLIST_ENV, raising=False)
    gitconfig = tmp_path / "gitconfig"
    gitconfig.write_text("")
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(gitconfig))
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    for var in (
        "EMAIL",
        "GIT_AUTHOR_NAME",
        "GIT_AUTHOR_EMAIL",
        "GIT_COMMITTER_NAME",
        "GIT_COMMITTER_EMAIL",
        "GIT_DIR",
        "GIT_INDEX_FILE",
        "GIT_WORK_TREE",
        "PRE_COMMIT_FROM_REF",
        "PRE_COMMIT_TO_REF",
    ):
        monkeypatch.delenv(var, raising=False)


@pytest.fixture
def repo(tmp_path, monkeypatch):
    path = tmp_path / "repo"
    path.mkdir()
    git(path, "init", "-q")
    monkeypatch.chdir(path)
    return path


# --- commit identities -------------------------------------------------------


@pytest.mark.parametrize(
    "email",
    [
        NOREPLY,
        "noreply@github.com",
        "49699333+dependabot[bot]@users.noreply.github.com",
        "noreply@anthropic.com",
    ],
)
def test_noreply_addresses_are_valid_commit_identities(email):
    assert guard.NOREPLY_EMAIL.match(email)


@pytest.mark.parametrize(
    "email",
    [
        addr("jane.doe", "mail.zz"),
        addr("x", "users.noreply.github.com.mail.zz"),
        "someone@example.org",
    ],
)
def test_other_addresses_are_not_valid_commit_identities(email):
    assert not guard.NOREPLY_EMAIL.match(email)


def test_identity_accepts_noreply_address(repo):
    git(repo, "config", "user.name", "Y3K5")
    git(repo, "config", "user.email", NOREPLY)
    assert guard.main(["identity"]) == 0


def test_identity_rejects_personal_address_without_printing_it(repo, capsys):
    email = addr("jane.doe", "mail.zz")
    git(repo, "config", "user.name", "Jane")
    git(repo, "config", "user.email", email)
    assert guard.main(["identity"]) == 1
    out = capsys.readouterr().out
    assert email not in out
    assert "j***@m***.zz" in out


def test_commits_flags_personal_address_without_printing_it(repo, capsys):
    commit(repo, NOREPLY)
    assert guard.main(["commits"]) == 0
    email = addr("jane.doe", "mail.zz")
    commit(repo, email)
    assert guard.main(["commits"]) == 1
    out = capsys.readouterr().out
    assert "author email" in out
    assert "committer email" in out
    assert email not in out


def test_commits_in_empty_repository_pass(repo):
    assert guard.main(["commits"]) == 0


def test_outgoing_commits_use_the_range_from_pre_commit(repo, monkeypatch):
    base = commit(repo, addr("jane.doe", "mail.zz"))  # pretend this one is already pushed
    head = commit(repo, NOREPLY)
    monkeypatch.setenv("PRE_COMMIT_FROM_REF", base)
    monkeypatch.setenv("PRE_COMMIT_TO_REF", head)
    assert guard.main(["commits", "--outgoing"]) == 0


def test_outgoing_commits_default_to_everything_not_on_a_remote(repo):
    commit(repo, addr("jane.doe", "mail.zz"))
    assert guard.main(["commits", "--outgoing"]) == 1


def test_denylisted_author_name_is_flagged(repo, monkeypatch, capsys):
    monkeypatch.setenv(guard.DENYLIST_ENV, SECRET_NAME)
    commit(repo, NOREPLY, name=SECRET_NAME)
    assert guard.main(["commits"]) == 1
    assert SECRET_NAME.casefold() not in capsys.readouterr().out.casefold()


# --- text contents -----------------------------------------------------------


def test_personal_email_is_flagged_and_masked():
    email = addr("jane.doe", "mail.zz")
    findings = guard.scan_text(f"Contact:\n{email}\n", "README.md", [])
    assert [finding.where for finding in findings] == ["README.md:2"]
    assert email not in str(findings[0])
    assert "j***@m***.zz" in str(findings[0])


@pytest.mark.parametrize(
    "text",
    [
        "someone@example.org",
        "git@github.com:Y3K5/MARSE.git",
        NOREPLY,
        "logo@2x.png",
        "rates = weights @ matrix.values",
    ],
)
def test_placeholder_addresses_are_allowed(text):
    assert guard.scan_text(text, "x", []) == []


@pytest.mark.parametrize(
    "text",
    [
        home("jane", "data.csv"),
        "file://" + home("jane"),
        "\\".join(["C:", "Users", "jane", "data.csv"]),
        "\\\\".join(["C:", "Users", "jane", "data.csv"]),  # JSON-escaped
        "/".join(["", "home", "jane", "runs"]),
    ],
)
def test_home_directory_paths_are_flagged(text):
    assert guard.scan_text(f'path = "{text}"', "x", [])


@pytest.mark.parametrize(
    "text",
    [
        "/home/runner/work/MARSE",
        "~/.config/marse/private-denylist.txt",
        "https://example.org/Users/jane",
        "/api/users/42",
        "C:\\Users\\Public\\data",
    ],
)
def test_generic_paths_are_allowed(text):
    assert guard.scan_text(text, "x", []) == []


def test_commit_message_checks_ignore_comments_and_the_diff(tmp_path):
    message = tmp_path / "COMMIT_EDITMSG"
    message.write_text(
        "Add growth model\n\n"
        "Co-Authored-By: Claude <noreply@anthropic.com>\n"
        "# Please enter the commit message for your changes.\n"
        "# ------------------------ >8 ------------------------\n"
        f"+ contact {addr('jane', 'mail.zz')}\n"
    )
    assert guard.main(["message", str(message)]) == 0
    message.write_text(f"Add growth model\n\nSigned-off-by: Jane <{addr('jane', 'mail.zz')}>\n")
    assert guard.main(["message", str(message)]) == 1


# --- files -------------------------------------------------------------------


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        (".env", "credential"),
        ("config/.env.local", "credential"),
        ("id_ed25519", "credential"),
        ("keys/server.pem", "credential"),
        ("MARSE_plan.pdf", "PDF"),
        ("notes.DOCX", "PDF"),
        ("data/reads.fastq.gz", "raw research data"),
        ("backup.zip", "archive"),
        (".DS_Store", "OS or log"),
        ("run.log", "OS or log"),
        ("private/notes.md", "private/"),
        ("Private/Plan.md", "private/"),
    ],
)
def test_sensitive_file_types_are_blocked(path, expected):
    reason = guard.blocked_reason(path)
    assert reason is not None
    assert expected in reason


@pytest.mark.parametrize(
    "path",
    [
        ".env.example",
        "README.md",
        "src/marse/core/__init__.py",
        "docs/figures/growth.png",
        "tests/data/golden.parquet",
    ],
)
def test_ordinary_files_are_not_blocked(path):
    assert guard.blocked_reason(path) is None


def test_files_command_reports_blocked_file(repo, capsys):
    (repo / "plan.pdf").write_bytes(b"%PDF-1.7")
    assert guard.main(["files", "plan.pdf"]) == 1
    assert "plan.pdf: never commit" in capsys.readouterr().out


def test_files_all_checks_tracked_files_only(repo):
    (repo / "notes.md").write_text("# Notes\n")
    (repo / "untracked.pdf").write_bytes(b"%PDF-1.7")
    git(repo, "add", "notes.md")
    assert guard.main(["files", "--all"]) == 0


def test_large_files_are_refused(repo, monkeypatch, capsys):
    monkeypatch.setattr(guard, "MAX_FILE_BYTES", 16)
    (repo / "trajectory.csv").write_text("t,biomass\n" * 4)
    assert guard.main(["files", "trajectory.csv"]) == 1
    assert "keep datasets" in capsys.readouterr().out


def _notebook(outputs: list) -> str:
    cell = {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": outputs,
        "source": "1 + 1",
    }
    return json.dumps({"cells": [cell], "metadata": {}, "nbformat": 4, "nbformat_minor": 5})


def test_notebooks_must_not_contain_outputs(repo):
    (repo / "demo.ipynb").write_text(_notebook([{"output_type": "stream", "text": "2"}]))
    (repo / "clean.ipynb").write_text(_notebook([]))
    assert guard.main(["files", "demo.ipynb"]) == 1
    assert guard.main(["files", "clean.ipynb"]) == 0


JPEG_WITH_EXIF = b"\xff\xd8\xff\xe1\x00\x16Exif\x00\x00MM\x00*\x00\x00\x00\x08\x00\x00\xff\xd9"
JPEG_WITHOUT_METADATA = (
    b"\xff\xd8\xff\xe0\x00\x10JFIF\x00\x01\x01\x00\x00\x01\x00\x01\x00\x00\xff\xd9"
)
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def test_photos_with_exif_metadata_are_flagged(repo):
    (repo / "photo.jpg").write_bytes(JPEG_WITH_EXIF)
    (repo / "plot.jpg").write_bytes(JPEG_WITHOUT_METADATA)
    assert guard.main(["files", "photo.jpg"]) == 1
    assert guard.main(["files", "plot.jpg"]) == 0


def test_xmp_is_flagged_only_when_it_holds_location_or_authorship():
    xmp = b'<x:xmpmeta xmlns:x="adobe:ns:meta/">%s</x:xmpmeta>'
    located = PNG_SIGNATURE + xmp % b"<exif:GPSLatitude>12,34N</exif:GPSLatitude>"
    screenshot = PNG_SIGNATURE + xmp % b"<exif:UserComment>Screenshot</exif:UserComment>"
    assert guard.image_metadata_problems(located, ".png")
    assert guard.image_metadata_problems(screenshot, ".png") == []


@pytest.mark.skipif(sys.platform == "win32", reason="creating symlinks needs extra privileges")
def test_symlink_targets_are_checked(repo):
    (repo / "datasets").symlink_to(home("jane", "datasets"))
    assert guard.main(["files", "datasets"]) == 1


# --- private denylist and allowlist -----------------------------------------


def test_denylist_matches_are_reported_but_never_echoed(repo, monkeypatch, capsys):
    monkeypatch.setenv(guard.DENYLIST_ENV, f"# my details\n{SECRET_NAME}\n")
    (repo / "README.md").write_text(f"Written by {SECRET_NAME.upper()}\n")
    (repo / f"{SECRET_NAME} notes.md").write_text("nothing personal\n")
    assert guard.main(["files", "README.md", f"{SECRET_NAME} notes.md"]) == 1
    out = capsys.readouterr().out
    assert "README.md:1: matches an entry in your private denylist" in out
    assert "*** notes.md" in out
    assert SECRET_NAME.casefold() not in out.casefold()


def test_denylist_file_is_loaded_from_the_config_directory(tmp_path, monkeypatch):
    monkeypatch.delenv(guard.DENYLIST_FILE_ENV)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    listing = tmp_path / "marse" / "private-denylist.txt"
    listing.parent.mkdir()
    listing.write_text(
        f"\ufeff# comments are ignored\n  {SECRET_NAME}  \n\n+1 555 0100\n", encoding="utf-8"
    )
    assert guard.load_denylist(os.environ) == ["+1 555 0100", SECRET_NAME.casefold()]


def test_allowlist_exempts_generic_checks_but_never_the_denylist(repo, monkeypatch):
    (repo / ".privacy-allowlist").write_text("# reviewed exceptions\ndocs/*.pdf\n")
    (repo / "docs").mkdir()
    poster = repo / "docs" / "poster.pdf"
    poster.write_bytes(b"%PDF-1.7 no metadata")
    assert guard.main(["files", "docs/poster.pdf"]) == 0
    monkeypatch.setenv(guard.DENYLIST_ENV, SECRET_NAME)
    poster.write_bytes(b"%PDF-1.7 /Author (" + SECRET_NAME.encode() + b")")
    assert guard.main(["files", "docs/poster.pdf"]) == 1
