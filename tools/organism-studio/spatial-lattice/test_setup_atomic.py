"""Invented structures qualify transport/topology, never biological interactions."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

import pytest
from setup_atomic import (
    build,
    build_topology,
    checked,
    digest,
    parse_pdb,
    qualify_atoms,
    qualify_pdb_prediction,
)


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
    def prediction_fixture(self):
        def fixed(fields):
            chars = [" "] * 80
            for start, width, value in fields:
                chars[start : start + width] = str(value).rjust(width)
            return "".join(chars)

        dbref = fixed(
            [
                (0, 6, "DBREF "),
                (7, 4, "XXXX"),
                (12, 1, "A"),
                (14, 4, 1),
                (20, 4, 2),
                (26, 6, "UNP"),
                (33, 8, "TEST1"),
                (42, 12, "TEST1_TEST"),
                (55, 5, 1),
                (62, 5, 2),
            ]
        )
        lines = [
            "TITLE     ALPHAFOLD MONOMER V2.0 PREDICTION FOR INVENTED (TEST1)",
            dbref,
            fixed([(0, 6, "SEQRES"), (7, 3, 1), (11, 1, "A"), (13, 4, 2)])[:19] + "ALA GLY",
        ]
        for n, name in enumerate(["ALA", "GLY"], 1):
            for i, (atom_name, element) in enumerate(
                [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]
            ):
                lines.append(
                    fixed(
                        [
                            (0, 6, "ATOM  "),
                            (6, 5, (n - 1) * 4 + i + 1),
                            (12, 4, atom_name),
                            (17, 3, name),
                            (21, 1, "A"),
                            (22, 4, n),
                            (30, 8, f"{n * 4 + i * 0.5:.3f}"),
                            (38, 8, "0.000"),
                            (46, 8, "0.000"),
                            (54, 6, "1.00"),
                            (60, 6, "80.00"),
                            (76, 2, element),
                        ]
                    )
                )
        raw = ("\n".join(lines) + "\n").encode()
        atoms, _ = parse_pdb(raw)
        shape = dict(
            id="AF-TEST1-F1",
            kind="prediction",
            admission="prediction_sequence_matched",
            chains=[dict(chain="A", ca=[[1, 4.5, 0, 0], [2, 8.5, 0, 0]])],
        )
        ledger = dict(chains=[dict(chain="A", residues=[[1, "ALA", ""], [2, "GLY", ""]])])
        protein = dict(
            accession="QUERY1",
            reference=dict(sequence="AG", length=2, sequence_sha256=digest(b"AG")),
            subtractive_record=dict(strain="invented", sequence_sha256=digest(b"AG")),
        )
        context = dict(id="AF-TEST1-F1", kind="prediction")
        return raw, atoms, shape, ledger, protein, context

    def test_prediction_requires_exact_sequence_model_and_original_numbering(self):
        raw, atoms, shape, ledger, protein, context = self.prediction_fixture()
        result = qualify_pdb_prediction(raw, atoms, shape, ledger, protein, context)
        assert result["sequence_sha256"] == digest(b"AG")
        assert atoms[0]["label_seq"] is None
        assert atoms[0]["charge"] is None
        assert atoms[0]["polymer_position"] == 1
        for change in [
            lambda p: p[4]["reference"].update(sequence="AA"),
            lambda p: p[4]["subtractive_record"].update(sequence_sha256="0" * 64),
            lambda p: p[5].update(kind="experimental"),
            lambda p: p[5].update(id="AF-HOMOLOG-F1"),
            lambda p: p[2].update(admission="related_prediction"),
            lambda p: p[1][0].update(insertion="A"),
            lambda p: p[1][1].update(xyz=[40, 0, 0]),
            lambda p: p[1].pop(),
        ]:
            values = copy.deepcopy(self.prediction_fixture())
            change(values)
            with pytest.raises(ValueError, match=r"[Pp]rediction|Atom/backbone"):
                qualify_pdb_prediction(*values)
        long = raw.replace(b"TEST1", b"A0A1TEST00")
        extended = copy.deepcopy(self.prediction_fixture())
        extended[2]["id"] = "AF-A0A1TEST00-F1"
        extended[5]["id"] = "AF-A0A1TEST00-F1"
        assert qualify_pdb_prediction(long, *extended[1:])["source_accession"] == "A0A1TEST00"
        for before, after in [
            (b"(TEST1)", b"(OTHER)"),
            (b"ALA GLY", b"ALA ALA"),
            (b"TEST1   ", b"WRONG   "),
        ]:
            changed = raw.replace(before, after)
            assert changed != raw
            with pytest.raises(ValueError, match=r"[Pp]rediction|Atom/backbone"):
                qualify_pdb_prediction(changed, *self.prediction_fixture()[1:])

    def test_prediction_import_preserves_inputs_and_records_missing_structures(self):
        raw, _, shape, entry, protein, context = self.prediction_fixture()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            runtime = root / "runtime"
            (runtime / "data").mkdir(parents=True)
            shape["source_view_sha256"] = digest(raw)
            protein["id"] = "invented"
            data = json.dumps(dict(proteins=[protein], shape_traces=dict(invented=shape))).encode()
            ledger = json.dumps(dict(entries=dict(invented=entry))).encode()
            for name, content in [("data/map.json", data), ("residues.json", ledger)]:
                (runtime / name).write_bytes(content)
            (runtime / "SOURCE_MANIFEST.json").write_text(json.dumps(dict(sha256=digest(data))))
            (runtime / "RESIDUE_MANIFEST.json").write_text(
                json.dumps(dict(sha256=digest(ledger), source_snapshot_sha256=digest(data)))
            )
            source = root / "prediction.pdb"
            source.write_bytes(raw)
            template = root / "template.cif"
            template.write_text("data_invented\n_entry.id INVENTED\n")
            record = dict(
                protein_id="invented",
                source_format="pdb_prediction",
                view_path=str(source),
                view_sha256=digest(raw),
                raw_path=str(source),
                raw_sha256=digest(raw),
                context=context,
            )
            config = dict(
                schema="marse-atomic-import-config-v1",
                runtime=str(runtime),
                output=str(root / "output"),
                template_path=str(template),
                template_sha256=digest(template.read_bytes()),
                models=[record],
                requested_proteins=["invented", "missing"],
            )
            result = build(config)
            assert set(result["entries"]) == {"invented"}
            assert result["holds"][0]["protein_id"] == "missing"
            payload = json.loads(
                (root / "output" / result["entries"]["invented"]["file"]).read_bytes()
            )
            assert len(payload["atoms"]) == 8
            assert payload["source_format"] == "pdb_prediction"
            assert payload["qualification"]["sequence_sha256"] == digest(b"AG")
            assert source.read_bytes() == raw
            assert (runtime / "data/map.json").read_bytes() == data
            assert (runtime / "residues.json").read_bytes() == ledger
            with pytest.raises(ValueError, match="already exists"):
                build(config)
            for mode in ["different_bytes", "unknown_format"]:
                config["output"] = str(root / mode)
                if mode == "different_bytes":
                    altered = root / "different.pdb"
                    altered.write_bytes(raw + b"\n")
                    record.update(raw_path=str(altered), raw_sha256=digest(altered.read_bytes()))
                else:
                    record.update(source_format="guessed")
                held = build(config)
                assert held["entries"] == {}
                assert held["holds"][0]["protein_id"] == "invented"
            assert source.read_bytes() == raw

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
