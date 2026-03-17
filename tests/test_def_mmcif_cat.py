"""Tests for pdbe_sifts.sifts_to_mmcif.def_mmcif_cat constants."""

from pdbe_sifts.sifts_to_mmcif.def_mmcif_cat import (
    DELTA_CAT,
    MAPPED_DB,
    NEW_MMCIF_CAT,
    PRIMARY_KEYS,
    SIFTS_ATOMSITE_ITEM,
    SIFTS_NEW_CAT,
)


class TestSiftsNewCat:
    def test_contains_atom_site(self):
        assert "atom_site" in SIFTS_NEW_CAT

    def test_contains_unp_segments(self):
        assert "pdbx_sifts_unp_segments" in SIFTS_NEW_CAT

    def test_is_list(self):
        assert isinstance(SIFTS_NEW_CAT, list)


class TestSiftsAtomsiteItem:
    def test_contains_id(self):
        assert "id" in SIFTS_ATOMSITE_ITEM

    def test_is_list(self):
        assert isinstance(SIFTS_ATOMSITE_ITEM, list)


class TestNewMmcifCat:
    def test_has_three_categories(self):
        assert len(NEW_MMCIF_CAT) == 3

    def test_unp_segments_key(self):
        assert "_pdbx_sifts_unp_segments" in NEW_MMCIF_CAT

    def test_unp_segments_has_unp_acc(self):
        assert "unp_acc" in NEW_MMCIF_CAT["_pdbx_sifts_unp_segments"]

    def test_xref_db_has_seq_id(self):
        assert "seq_id" in NEW_MMCIF_CAT["_pdbx_sifts_xref_db"]

    def test_all_values_are_lists(self):
        for v in NEW_MMCIF_CAT.values():
            assert isinstance(v, list)
            assert len(v) > 0


class TestPrimaryKeys:
    def test_atom_site_has_id(self):
        assert PRIMARY_KEYS["_atom_site"] == ["id"]

    def test_unp_segments_has_unp_acc(self):
        assert "unp_acc" in PRIMARY_KEYS["_pdbx_sifts_unp_segments"]

    def test_all_values_non_empty(self):
        for v in PRIMARY_KEYS.values():
            assert len(v) > 0


class TestDeltaCat:
    def test_has_two_categories(self):
        assert len(DELTA_CAT) == 2

    def test_unp_segments_has_ranges(self):
        cols = DELTA_CAT["_pdbx_sifts_unp_segments"]
        assert "unp_start" in cols
        assert "unp_end" in cols


class TestMappedDb:
    def test_contains_unp(self):
        assert "UNP" in MAPPED_DB

    def test_contains_pfam(self):
        assert "Pfam" in MAPPED_DB

    def test_is_list(self):
        assert isinstance(MAPPED_DB, list)
