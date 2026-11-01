"""Invented transport checks; no campaign data or scientific run."""

import json
from pathlib import Path

import pytest
from fixture_relations import make_fixture, raw, sha
from setup_relations import install, read_build


def test_verified_json_only_and_idempotent_import(tmp_path: Path):
    source = tmp_path / "build"
    make_fixture(source)
    output = tmp_path / "runtime"
    selected = install(source, output)
    assert set(selected) == {
        "SOURCE_MANIFEST.json",
        "RESIDUE_MANIFEST.json",
        "residues.json",
        "data/map.json",
        "protein-relation-map.json",
        "RELATION_BUILD.json",
        "BUILD_MANIFEST.json",
        "IMPORT.json",
    }
    assert install(source / "studio", output) == selected
    assert not any(p.suffix == ".js" for p in output.rglob("*"))
    assert "source_path" not in json.loads((output / "IMPORT.json").read_text())


@pytest.mark.parametrize(
    "kind",
    ["changed", "missing", "traversal", "symlink", "scientific", "trace", "model", "binding"],
)
def test_bad_retained_build_refused_before_output(tmp_path: Path, kind: str):
    source = tmp_path / "build"
    manifest = make_fixture(source)
    path = source / "studio/data/map.json"
    if kind == "changed":
        path.write_bytes(b"changed")
    elif kind == "missing":
        path.unlink()
    elif kind == "traversal":
        manifest["files"]["../outside.json"] = "a" * 64
    elif kind == "symlink":
        outside = tmp_path / "outside.json"
        outside.write_bytes(path.read_bytes())
        path.unlink()
        path.symlink_to(outside)
    elif kind == "scientific":
        manifest["scientific_execution"] = True
    elif kind == "trace":
        manifest["study_revision_sha256"] = "a" * 64
    elif kind == "model":
        manifest["model_hashes"]["review-a"] = "a" * 64
    else:
        data = json.loads(path.read_text())
        data["shape_traces"]["review-a"]["source_view_sha256"] = "a" * 64
        path.write_bytes(raw(data))
        manifest["files"]["studio/data/map.json"] = sha(path.read_bytes())
    (source / "BUILD_MANIFEST.json").write_bytes(raw(manifest))
    output = tmp_path / "runtime"
    with pytest.raises((ValueError, OSError)):
        install(source, output)
    assert not output.exists()


def test_different_existing_runtime_is_preserved(tmp_path: Path):
    source = tmp_path / "build"
    make_fixture(source)
    output = tmp_path / "runtime"
    output.mkdir()
    private = output / "residues.json"
    private.write_bytes(b"preserved input")
    with pytest.raises(ValueError, match="Existing runtime differs"):
        install(source, output)
    assert private.read_bytes() == b"preserved input"


def test_extra_supplied_code_is_checked_but_never_installed(tmp_path: Path):
    source = tmp_path / "build"
    manifest = make_fixture(source)
    code = source / "studio" / "untrusted.js"
    code.write_bytes(b"throw Error('never execute');")
    manifest["files"]["studio/untrusted.js"] = sha(code.read_bytes())
    (source / "BUILD_MANIFEST.json").write_bytes(raw(manifest))
    selected = read_build(source)
    assert all(not name.endswith(".js") for name in selected)


def test_distribution_excludes_private_runtime_and_existing_asset_bundles():
    import tomllib

    root = Path(__file__).resolve().parents[3]
    data = tomllib.loads((root / "pyproject.toml").read_text())
    excluded = data["tool"]["hatch"]["build"]["targets"]["sdist"]["exclude"]
    assert "/tools/organism-studio/spatial-lattice/runtime/**" in excluded
    assert "/tools/organism-studio/spatial-lattice/.relations-import-*/**" in excluded
    assert "/tools/organism-studio/vendor/**" in excluded


def test_import_retains_the_bytes_that_passed_checksum(tmp_path: Path, monkeypatch):
    source = tmp_path / "build"
    make_fixture(source)
    payload = source / "studio/residues.json"
    verified = payload.read_bytes()
    original_read = Path.read_bytes
    reads = 0

    def changing_read(path):
        nonlocal reads
        if path == payload:
            reads += 1
            if reads > 1:
                return b"unchecked replacement"
        return original_read(path)

    monkeypatch.setattr(Path, "read_bytes", changing_read)
    selected = read_build(source)
    assert selected["residues.json"] == verified
    assert reads == 1
