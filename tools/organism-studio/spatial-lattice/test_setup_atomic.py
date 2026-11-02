"""Invented structures qualify transport/topology, never biological interactions."""

import json
import tempfile
import unittest
from pathlib import Path

import pytest
from setup_atomic import build, build_topology, checked, digest, parse_pdb, qualify_atoms


def column_doc(category, rows):
    return {category + "." + key: [row.get(key, "?") for row in rows] for key in rows[0]}


def atom(i, name, x, seq=1, residue=8):
    return dict(
        i=i,
        name=name,
        xyz=[x, 0, 0],
        chain="A",
        residue=residue,
        insertion="A" if seq == 1 else "",
        resname="ALA",
        element="C",
        group="ATOM",
        label_seq=seq,
        label_chain="A",
        altloc="",
    )


class AtomicChecks(unittest.TestCase):
    def test_no_proximity_bond_and_no_missing_atom_creation(self):
        atoms = [atom(0, "CA", 0), atom(1, "C", 1.4)]
        bonds, limits = build_topology(atoms, {}, {}, "raw", "fallback")
        assert bonds == []
        assert limits["holds"]
        doc = column_doc(
            "_chem_comp_bond",
            [
                dict(comp_id="ALA", atom_id_1="CA", atom_id_2="C", value_order="SING"),
                dict(comp_id="ALA", atom_id_1="CA", atom_id_2="H", value_order="SING"),
            ],
        )
        bonds, _ = build_topology(atoms, doc, {}, "raw", "fallback")
        assert len(atoms) == 2
        assert [(b["kind"], b["order"]) for b in bonds] == [("component_template", "sing")]
        atoms[1]["xyz"] = [10, 0, 0]
        bonds, limits = build_topology(atoms, doc, {}, "raw", "fallback")
        assert not bonds
        assert limits["holds"]

    def test_polymer_requires_declared_consecutive_positions_chain_and_geometry(self):
        atoms = [atom(0, "C", 0), atom(1, "N", 1.3, 2, 9)]
        bonds, _ = build_topology(atoms, {}, {}, "raw", "fallback")
        assert bonds[0]["kind"] == "polymer_link"
        for key, value in [("label_seq", 3), ("label_chain", "B"), ("xyz", [5, 0, 0])]:
            changed = [dict(a) for a in atoms]
            changed[1][key] = value
            assert not build_topology(changed, {}, {}, "raw", "fallback")[0]

    def test_explicit_connections_keep_type_and_symmetry_hold(self):
        atoms = [atom(0, "SG", 0), atom(1, "SG", 2, 2, 9)]
        row = dict(
            id="invented",
            conn_type_id="disulf",
            ptnr1_auth_asym_id="A",
            ptnr2_auth_asym_id="A",
            ptnr1_auth_seq_id="8",
            ptnr2_auth_seq_id="9",
            ptnr1_label_atom_id="SG",
            ptnr2_label_atom_id="SG",
            pdbx_ptnr1_PDB_ins_code="A",
            pdbx_ptnr2_PDB_ins_code="?",
            ptnr1_symmetry="1_555",
            ptnr2_symmetry="1_555",
        )
        for kind in ("disulf", "metalc", "hydrog", "covale"):
            row["conn_type_id"] = kind
            bonds, _ = build_topology(
                atoms, column_doc("_struct_conn", [row]), {}, "raw", "fallback"
            )
            assert bonds[0]["kind"] == "source_" + kind
        for key, value in [
            ("ptnr2_symmetry", "2_555"),
            ("ptnr2_auth_asym_id", "B"),
            ("pdbx_ptnr1_PDB_ins_code", "?"),
        ]:
            changed = dict(row, **{key: value})
            bonds, limits = build_topology(
                atoms, column_doc("_struct_conn", [changed]), {}, "raw", "fallback"
            )
            assert not bonds
            assert limits["holds"]

    def test_atom_source_and_backbone_join_rejects_wrong_identity(self):
        atoms = [atom(0, "CA", 0)]
        row = dict(
            auth_asym_id="A",
            auth_seq_id="8",
            pdbx_PDB_ins_code="A",
            auth_comp_id="ALA",
            auth_atom_id="CA",
            label_alt_id="?",
            Cartn_x="0",
            Cartn_y="0",
            Cartn_z="0",
            type_symbol="C",
            label_seq_id="1",
            label_asym_id="A",
        )
        shape = dict(chains=[dict(chain="A", ca=[[8, 0, 0, 0]])])
        ledger = dict(chains=[dict(chain="A", residues=[[8, "ALA", "A"]])])
        qualify_atoms(atoms, column_doc("_atom_site", [row]), shape, ledger)
        for key, value in [
            ("auth_seq_id", "32"),
            ("pdbx_PDB_ins_code", "?"),
            ("auth_asym_id", "B"),
            ("type_symbol", "N"),
            ("Cartn_x", ".1"),
        ]:
            with pytest.raises(ValueError, match="View/source atom"):
                qualify_atoms(
                    atoms, column_doc("_atom_site", [dict(row, **{key: value})]), shape, ledger
                )
        with pytest.raises(ValueError, match="View/source atom identity"):
            qualify_atoms(atoms, column_doc("_atom_site", [row, row]), shape, ledger)

    def test_water_context_reports_retained_atoms_without_reconstruction(self):
        atoms = [dict(atom(0, "O", 0), resname="HOH", group="HETATM", label_seq=None)]
        before = [dict(a) for a in atoms]
        bonds, limits = build_topology(atoms, {}, {}, "raw", "fallback")
        assert not bonds
        assert limits["water_atoms_present"] == 1
        assert "1 water atoms retained" in limits["water_scope"]
        assert "omits water" not in limits["water_scope"]
        assert len(atoms) == len(before)

    def test_fresh_import_pins_entry_identity_and_refuses_stale_runtime(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            runtime = root / "runtime"
            (runtime / "data").mkdir(parents=True)
            chars = [" "] * 80
            for start, width, value in [
                (0, 6, "ATOM  "),
                (6, 5, "1"),
                (12, 4, "CA"),
                (17, 3, "ALA"),
                (21, 1, "A"),
                (22, 4, "8"),
                (26, 1, "A"),
                (30, 8, "0.000"),
                (38, 8, "0.000"),
                (46, 8, "0.000"),
                (54, 6, "1.00"),
                (60, 6, "20.00"),
                (76, 2, "C"),
            ]:
                chars[start : start + width] = value.rjust(width)
            pdb = "".join(chars).encode()
            view = root / "view.pdb"
            view.write_bytes(pdb)
            columns = [
                "auth_asym_id",
                "auth_seq_id",
                "pdbx_PDB_ins_code",
                "auth_comp_id",
                "auth_atom_id",
                "label_alt_id",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
                "type_symbol",
                "label_seq_id",
                "label_asym_id",
            ]

            def cif(entry):
                return (
                    "data_invented\n_entry.id "
                    + entry
                    + "\nloop_\n"
                    + "\n".join("_atom_site." + c for c in columns)
                    + "\nA 8 A ALA CA . 0 0 0 C 1 A\n"
                ).encode()

            raw = root / "raw.cif"
            raw.write_bytes(cif("TEST"))
            shape = dict(
                source_view_sha256=digest(pdb), chains=[dict(chain="A", ca=[[8, 0, 0, 0]])]
            )
            data = json.dumps(dict(shape_traces=dict(invented=shape))).encode()
            ledger = json.dumps(
                dict(
                    entries=dict(
                        invented=dict(chains=[dict(chain="A", residues=[[8, "ALA", "A"]])])
                    )
                )
            ).encode()
            (runtime / "data/map.json").write_bytes(data)
            (runtime / "residues.json").write_bytes(ledger)
            (runtime / "SOURCE_MANIFEST.json").write_text(json.dumps(dict(sha256=digest(data))))
            ledger_meta = dict(sha256=digest(ledger), source_snapshot_sha256=digest(data))
            (runtime / "RESIDUE_MANIFEST.json").write_text(json.dumps(ledger_meta))
            template = root / "template.cif"
            template.write_bytes(cif("TEST"))
            record = dict(
                protein_id="invented",
                view_path=str(view),
                view_sha256=digest(pdb),
                raw_path=str(raw),
                raw_sha256=digest(raw.read_bytes()),
                context=dict(
                    id="TEST", kind="experimental", source_url="https://example.org/structure"
                ),
            )
            config = dict(
                schema="marse-atomic-import-config-v1",
                runtime=str(runtime),
                output=str(root / "output"),
                template_path=str(template),
                template_sha256=digest(template.read_bytes()),
                models=[record],
            )
            result = build(config)
            assert not result["holds"]
            entry = result["entries"]["invented"]
            assert entry["atoms"] == 1
            assert (root / "output" / entry["pdb_file"]).read_bytes() == pdb
            with pytest.raises(ValueError, match="already exists"):
                build(config)
            for n, source_id, context_id, kind, admitted in [
                (1, "TEST-F1", "TEST", "experimental", False),
                (2, "TEST-F1", "TEST", "prediction", False),
                (3, "AF-INVENTED-F1", "AF-INVENTED", "prediction", True),
            ]:
                raw.write_bytes(cif(source_id))
                record["raw_sha256"] = digest(raw.read_bytes())
                record["context"].update(id=context_id, kind=kind)
                config["output"] = str(root / ("output" + str(n)))
                result = build(config)
                assert ("invented" in result["entries"]) == admitted
                assert bool(result["holds"]) != admitted
            ledger_meta["source_snapshot_sha256"] = "0" * 64
            (runtime / "RESIDUE_MANIFEST.json").write_text(json.dumps(ledger_meta))
            config["output"] = str(root / "stale")
            with pytest.raises(ValueError, match="source/residue pins differ"):
                build(config)
            assert not (root / "stale").exists()
            assert view.read_bytes() == pdb

    def test_hash_and_immutable_destination(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "source"
            path.write_bytes(b"invented")
            assert checked(str(path), digest(b"invented")) == b"invented"
            with pytest.raises(ValueError, match="checksum or size differs"):
                checked(str(path), "0" * 64)
            with pytest.raises(ValueError, match="already exists"):
                build(dict(schema="marse-atomic-import-config-v1", runtime=tmp, output=tmp))
            assert path.read_bytes() == b"invented"
            dangling = Path(tmp) / "dangling"
            dangling.symlink_to(Path(tmp) / "missing")
            with pytest.raises(ValueError, match="already exists"):
                build(
                    dict(schema="marse-atomic-import-config-v1", runtime=tmp, output=str(dangling))
                )
            assert dangling.is_symlink()
        with pytest.raises(ValueError, match="No retained primary atom"):
            parse_pdb(b"END\n")


if __name__ == "__main__":
    unittest.main()
