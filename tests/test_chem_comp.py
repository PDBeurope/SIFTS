"""Tests for pdbe_sifts.mmcif.chem_comp.ChemCompMapping."""

import pytest
from pdbe_sifts.mmcif.chem_comp import STANDARD_AA, ChemCompMapping


@pytest.fixture(scope="module")
def cc():
    return ChemCompMapping()


class TestStandardAA:
    def test_has_20_entries(self):
        assert len(STANDARD_AA) == 20

    def test_known_mappings(self):
        assert STANDARD_AA["ALA"] == "A"
        assert STANDARD_AA["GLY"] == "G"
        assert STANDARD_AA["TRP"] == "W"
        assert STANDARD_AA["MET"] == "M"


class TestChemCompMappingGet:
    def test_standard_aa(self, cc):
        assert cc.get("ALA") == "A"
        assert cc.get("GLY") == "G"
        assert cc.get("TRP") == "W"

    def test_case_insensitive(self, cc):
        assert cc.get("ala") == "A"
        assert cc.get("Gly") == "G"

    def test_unknown_returns_X(self, cc):
        assert cc.get("XYZ") == "X"

    def test_list_input(self, cc):
        result = cc.get(["ALA", "GLY", "TRP"])
        assert result == "AGW"

    def test_none_input_returns_none(self, cc):
        assert cc.get(None) is None


class TestChemCompMappingCheck:
    def test_all_standard_aa_present(self, cc):
        # _check_chem_comp should not raise
        cc._check_chem_comp()
