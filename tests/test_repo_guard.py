"""Tests for tools/repo_guard.py.

A guard that never fires is worse than none, because it looks like protection.
Each rule is tested against a workflow that violates it and one that does not.
"""

from __future__ import annotations

import importlib.util
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PINNED = "actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1"


def _load_guard():
    spec = importlib.util.spec_from_file_location("repo_guard", ROOT / "tools" / "repo_guard.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


guard = _load_guard()


def workflow(body: str, tmp_path: Path, name: str = "w.yml") -> Path:
    path = tmp_path / name
    path.write_text(textwrap.dedent(body).lstrip(), encoding="utf-8")
    return path


def problems(body: str, tmp_path: Path) -> list[str]:
    return [p.detail for p in guard.check_workflow(workflow(body, tmp_path), "w.yml")]


SOUND = """
    name: CI
    on: [push]
    permissions:
      contents: read
    jobs:
      build:
        runs-on: ubuntu-latest
        steps:
          - uses: {pinned}
          - run: python -m pytest
""".replace("{pinned}", PINNED)


def test_the_repository_itself_passes_every_rule():
    assert guard.main([str(ROOT)]) == 0


def test_a_sound_workflow_raises_nothing(tmp_path):
    assert problems(SOUND, tmp_path) == []


# --- pinning -----------------------------------------------------------------


@pytest.mark.parametrize(
    "ref",
    ["actions/checkout@v4", "actions/checkout@main", "actions/checkout", "actions/checkout@abc123"],
)
def test_unpinned_actions_are_rejected(ref, tmp_path):
    found = problems(SOUND.replace(PINNED, ref), tmp_path)
    assert any("not pinned" in p for p in found), found


def test_local_and_docker_actions_are_exempt(tmp_path):
    for ref in ("./.github/actions/setup", "docker://alpine:3.20"):
        assert not any("not pinned" in p for p in problems(SOUND.replace(PINNED, ref), tmp_path))


# --- dangerous triggers ------------------------------------------------------


def test_pull_request_target_is_rejected(tmp_path):
    found = problems(SOUND.replace("on: [push]", "on:\n  pull_request_target:"), tmp_path)
    assert any("pull_request_target" in p for p in found), found


def test_ordinary_pull_request_is_allowed(tmp_path):
    found = problems(SOUND.replace("on: [push]", "on:\n  pull_request:"), tmp_path)
    assert found == []


# --- permissions -------------------------------------------------------------


def test_write_permissions_are_rejected_unless_declared(tmp_path):
    found = problems(SOUND.replace("contents: read", "contents: write"), tmp_path)
    assert any("contents: write" in p for p in found), found


def test_a_declared_write_permission_is_allowed(tmp_path, monkeypatch):
    monkeypatch.setitem(guard.ALLOWED_WRITE_PERMISSIONS, "w.yml", {"contents"})
    assert problems(SOUND.replace("contents: read", "contents: write"), tmp_path) == []


def test_a_workflow_without_a_permissions_block_is_rejected(tmp_path):
    body = SOUND.replace("    permissions:\n      contents: read\n", "")
    assert "permissions" not in body, "the fixture must really have lost its permissions block"
    assert any("least privilege" in p for p in problems(body, tmp_path))


def test_a_commented_out_write_permission_is_not_flagged(tmp_path):
    body = SOUND.replace("contents: read", "contents: read # not contents: write")
    assert problems(body, tmp_path) == []


# --- script injection --------------------------------------------------------


def test_untrusted_interpolation_into_a_shell_command_is_rejected(tmp_path):
    body = SOUND.replace(
        "- run: python -m pytest",
        '- run: echo "${{ github.event.pull_request.title }}"',
    )
    found = problems(body, tmp_path)
    assert any("shell command" in p for p in found), found


def test_untrusted_interpolation_inside_a_block_scalar_is_rejected(tmp_path):
    body = SOUND.replace(
        "- run: python -m pytest",
        "- run: |\n              echo ${{ github.head_ref }}\n              make build",
    )
    assert any("shell command" in p for p in problems(body, tmp_path))


def test_an_expression_in_a_with_input_is_not_a_shell_injection(tmp_path):
    """The real privacy.yml pattern: a ref passed to checkout, not to a shell."""
    body = """
        name: CI
        on: [push]
        permissions:
          contents: read
        jobs:
          build:
            runs-on: ubuntu-latest
            steps:
              - uses: PINNED
                with:
                  ref: ${{ github.event.pull_request.head.sha || github.sha }}
              - run: python3 tools/privacy_guard.py commits
              - uses: PINNED
                with:
                  ref: ${{ github.event.pull_request.head.sha || github.sha }}
              - run: python3 tools/privacy_guard.py files --all
    """.replace("PINNED", PINNED)
    assert problems(body, tmp_path) == []


def test_trusted_context_in_a_shell_command_is_allowed(tmp_path):
    body = SOUND.replace("- run: python -m pytest", "- run: echo ${{ github.sha }}")
    assert problems(body, tmp_path) == []


def test_an_env_indirection_is_the_accepted_workaround(tmp_path):
    body = SOUND.replace(
        "- run: python -m pytest",
        "- env:\n              TITLE: ${{ github.event.pull_request.title }}\n"
        '            run: echo "$TITLE"',
    )
    assert problems(body, tmp_path) == []


# --- publishing --------------------------------------------------------------


@pytest.mark.parametrize(
    "step",
    [
        "- uses: pypa/gh-action-pypi-publish@0ab0b79471669eb3a4d647e625009c62f9f3b241",
        "- run: twine upload dist/*",
        "- run: gh release create v1.0.0",
    ],
)
def test_publishing_without_an_environment_gate_is_rejected(step, tmp_path):
    body = SOUND.replace("- run: python -m pytest", step)
    found = problems(body, tmp_path)
    assert any("environment" in p for p in found), found


def test_publishing_behind_an_environment_is_allowed(tmp_path):
    body = SOUND.replace(
        "    runs-on: ubuntu-latest",
        "    runs-on: ubuntu-latest\n    environment: release",
    ).replace("- run: python -m pytest", "- run: twine upload dist/*")
    assert problems(body, tmp_path) == []


# --- structure ---------------------------------------------------------------


def test_a_missing_required_file_is_reported(tmp_path):
    (tmp_path / "src" / "marse").mkdir(parents=True)
    found = [str(p) for p in guard.check_structure(tmp_path, [])]
    assert any("README.md: required file is missing" in p for p in found), found


def test_a_removed_architecture_subpackage_is_reported(tmp_path):
    package = tmp_path / "src" / "marse"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text('"""Package."""\n')
    found = [str(p) for p in guard.check_structure(tmp_path, [])]
    assert any("src/marse/core" in p for p in found), found


def test_a_module_without_a_docstring_is_reported(tmp_path):
    package = tmp_path / "src" / "marse"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("X = 1\n")
    found = [str(p) for p in guard.check_structure(tmp_path, [])]
    assert any("no docstring" in p for p in found), found


def test_sprawl_at_the_repository_root_is_reported(tmp_path):
    found = [str(p) for p in guard.check_structure(tmp_path, ["scratch_notes.py", "README.md"])]
    assert any("scratch_notes.py" in p and "unexpected entry" in p for p in found), found


def test_known_root_entries_are_not_reported(tmp_path):
    found = [
        str(p) for p in guard.check_structure(tmp_path, ["src/marse/__init__.py", "docs/x.md"])
    ]
    assert not any("unexpected entry" in p for p in found)


def test_a_repository_with_no_workflows_is_reported(tmp_path):
    found = [str(p) for p in guard.check_workflows(tmp_path)]
    assert any("no workflows" in p for p in found), found


def test_main_reports_problems_and_exits_non_zero(tmp_path, capsys):
    assert guard.main([str(tmp_path)]) == 1
    assert "GOVERNANCE.md" in capsys.readouterr().out
