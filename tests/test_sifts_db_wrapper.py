"""Tests for pdbe_sifts.database.sifts_db_wrapper."""

from types import SimpleNamespace

import duckdb
import pytest
from pdbe_sifts.database.sifts_db_wrapper import (
    SiftsDB,
    residue_to_dict,
    segment_to_dict,
)


def _make_segment(**kwargs):
    defaults = {
        "entry_id": "1cbs",
        "entity_id": 1,
        "id": 1,
        "auth_asym_id": "A",
        "struct_asym_id": "A",
        "accession": "P29373",
        "name": "CRABP-II",
        "seq_version": 1,
        "unp_start": 1,
        "pdb_start": 1,
        "unp_end": 137,
        "pdb_end": 137,
        "auth_start": 1,
        "auth_start_icode": None,
        "auth_end": 137,
        "auth_end_icode": None,
        "conflicts": 0,
        "modifications": None,
        "unp_alignment": "AAAAAA",
        "pdb_alignment": "AAAAAA",
        "identity": 1.0,
        "score": 100.0,
        "best_mapping": True,
        "canonical_acc": True,
        "reference_acc": None,
        "chimera": False,
    }
    defaults.update(kwargs)
    return SimpleNamespace(**defaults)


def _make_residue(**kwargs):
    defaults = {
        "entry_id": "1cbs",
        "entity_id": 1,
        "id": 1,
        "auth_asym_id": "A",
        "struct_asym_id": "A",
        "unp_segment_id": 1,
        "auth_seq_id": 1,
        "auth_seq_id_ins_code": None,
        "pdb_seq_id": 1,
        "unp_seq_id": 1,
        "observed": "Y",
        "dbentry_id": 12345,
        "accession": "P29373",
        "name": "CRABP-II",
        "type": "ATOM",
        "unp_one_letter_code": "A",
        "pdb_one_letter_code": "A",
        "chem_comp_id": "ALA",
        "mh_id": None,
        "tax_id": 9606,
        "canonical_acc": True,
        "reference_acc": None,
        "best_mapping": True,
        "residue_id": "1cbs_A_1",
    }
    defaults.update(kwargs)
    return SimpleNamespace(**defaults)


class TestSegmentToDict:
    def test_returns_dict(self):
        seg = _make_segment()
        d = segment_to_dict(seg)
        assert isinstance(d, dict)

    def test_entry_id(self):
        seg = _make_segment(entry_id="1abc")
        assert segment_to_dict(seg)["entry_id"] == "1abc"

    def test_entity_id_is_int(self):
        seg = _make_segment(entity_id="2")
        d = segment_to_dict(seg)
        assert isinstance(d["entity_id"], int)
        assert d["entity_id"] == 2

    def test_best_mapping_is_bool(self):
        seg = _make_segment(best_mapping=1)
        d = segment_to_dict(seg)
        assert isinstance(d["best_mapping"], bool)

    def test_chimera_is_bool(self):
        seg = _make_segment(chimera=0)
        d = segment_to_dict(seg)
        assert d["chimera"] is False

    def test_all_expected_keys_present(self):
        expected = {
            "entry_id",
            "entity_id",
            "id",
            "auth_asym_id",
            "struct_asym_id",
            "accession",
            "identity",
            "score",
            "best_mapping",
            "canonical_acc",
            "chimera",
        }
        d = segment_to_dict(_make_segment())
        assert expected.issubset(d.keys())


class TestResidueToDict:
    def test_returns_dict(self):
        r = _make_residue()
        d = residue_to_dict(r)
        assert isinstance(d, dict)

    def test_entry_id(self):
        r = _make_residue(entry_id="2xyz")
        assert residue_to_dict(r)["entry_id"] == "2xyz"

    def test_entity_id_is_int(self):
        r = _make_residue(entity_id="3")
        d = residue_to_dict(r)
        assert isinstance(d["entity_id"], int)

    def test_best_mapping_none_preserved(self):
        r = _make_residue(best_mapping=None)
        d = residue_to_dict(r)
        assert d["best_mapping"] is None

    def test_best_mapping_bool_coercion(self):
        r = _make_residue(best_mapping=1)
        d = residue_to_dict(r)
        assert d["best_mapping"] is True

    def test_canonical_acc_is_bool(self):
        r = _make_residue(canonical_acc=0)
        d = residue_to_dict(r)
        assert d["canonical_acc"] is False


class TestSiftsDB:
    @pytest.fixture
    def db(self):
        conn = duckdb.connect(":memory:")
        return SiftsDB(conn), conn

    def test_tables_created(self, db):
        sdb, conn = db
        tables = {row[0] for row in conn.execute("SHOW TABLES").fetchall()}
        assert "sifts_xref_segment" in tables
        assert "sifts_xref_residue" in tables
