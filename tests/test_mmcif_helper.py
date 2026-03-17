"""Tests for pdbe_sifts.mmcif.mmcif_helper using the 1cbs.cif fixture."""

from pathlib import Path

import pytest
from pdbe_sifts.mmcif.chem_comp import ChemCompMapping
from pdbe_sifts.mmcif.mmcif_helper import mmCIF

CIF_1CBS = Path(__file__).parent / "data" / "cif" / "1cbs.cif"


@pytest.fixture(scope="module")
def mmcif_1cbs():
    cc = ChemCompMapping()
    return mmCIF("1cbs", cc, str(CIF_1CBS))


class TestMmcifConstruction:
    def test_loads_without_error(self, mmcif_1cbs):
        assert mmcif_1cbs is not None

    def test_pdbid(self, mmcif_1cbs):
        assert mmcif_1cbs.pdbid == "1cbs"

    def test_poly_seq_loaded(self, mmcif_1cbs):
        assert mmcif_1cbs.poly_seq is not None

    def test_file_not_found_raises(self):
        cc = ChemCompMapping()
        with pytest.raises(FileNotFoundError):
            mmCIF("XXXX", cc, "/nonexistent/path/XXXX.cif")


class TestGetChains:
    def test_returns_dict(self, mmcif_1cbs):
        chains = mmcif_1cbs.get_chains()
        assert isinstance(chains, dict)

    def test_chain_A_present(self, mmcif_1cbs):
        chains = mmcif_1cbs.get_chains()
        assert "A" in chains

    def test_chain_has_entity_and_asym(self, mmcif_1cbs):
        chains = mmcif_1cbs.get_chains()
        for entity_id, asym_id in chains.values():
            assert entity_id is not None
            assert asym_id is not None


class TestGetSequence:
    def test_returns_string(self, mmcif_1cbs):
        chains = mmcif_1cbs.get_chains()
        entity_id = next(iter(chains.values()))[0]
        seq = mmcif_1cbs.get_sequence(entity_id)
        assert isinstance(seq, str)

    def test_sequence_not_empty(self, mmcif_1cbs):
        chains = mmcif_1cbs.get_chains()
        entity_id = next(iter(chains.values()))[0]
        seq = mmcif_1cbs.get_sequence(entity_id)
        assert len(seq) > 0

    def test_sequence_is_uppercase_letters(self, mmcif_1cbs):
        chains = mmcif_1cbs.get_chains()
        entity_id = next(iter(chains.values()))[0]
        seq = mmcif_1cbs.get_sequence(entity_id)
        assert all(c.isupper() or c == "X" for c in seq)


class TestGetEntityId:
    def test_returns_entity_id_for_chain_A(self, mmcif_1cbs):
        eid = mmcif_1cbs.get_entity_id("A")
        assert eid is not None

    def test_unknown_chain_returns_none(self, mmcif_1cbs):
        eid = mmcif_1cbs.get_entity_id("Z")
        assert eid is None


class TestGetUnp:
    def test_returns_list(self, mmcif_1cbs):
        unps = mmcif_1cbs.get_unp("A")
        assert isinstance(unps, list)


class TestIsPoly:
    def test_entity_1_is_polypeptide(self, mmcif_1cbs):
        assert mmcif_1cbs.is_poly("1") is True


class TestGetResidues:
    def test_returns_non_empty_list(self, mmcif_1cbs):
        residues = mmcif_1cbs.get_residues("A")
        assert len(residues) > 0

    def test_each_item_is_int_tuple(self, mmcif_1cbs):
        residues = mmcif_1cbs.get_residues("A")
        for n, _ in residues:
            assert isinstance(n, int)
