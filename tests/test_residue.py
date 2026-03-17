"""Tests for pdbe_sifts.mmcif.residue.Residue."""

from pdbe_sifts.mmcif.residue import Residue


def _make_residue(**kwargs):
    defaults = {
        "n": 1,
        "auth_n": 1,
        "auth_ins": False,
        "oneL": "A",
        "threeL": "ALA",
        "rtype": "ATOM",
        "observed": True,
        "oneL_original": "A",
        "threeL_original": "ALA",
        "mh": None,
    }
    defaults.update(kwargs)
    return Residue(**defaults)


class TestResidueConstruction:
    def test_basic_attributes(self):
        r = _make_residue()
        assert r.n == 1
        assert r.auth_n == 1
        assert r.auth_ins is False
        assert r.oneL == "A"
        assert r.threeL == "ALA"
        assert r.rtype == "ATOM"
        assert r.observed is True
        assert r.mh is None

    def test_is_chromophore_false(self):
        r = _make_residue(rtype="ATOM")
        assert r.is_chromophore is False

    def test_is_chromophore_true(self):
        r = _make_residue(rtype="Chromophore")
        assert r.is_chromophore is True

    def test_oneL_original_and_threeL_original(self):
        r = _make_residue(oneL_original="G", threeL_original="GLY")
        assert r.oneL_original == "G"
        assert r.threeL_original == "GLY"


class TestResidueRepr:
    def test_repr_contains_n(self):
        r = _make_residue(n=42)
        assert "42" in repr(r)

    def test_repr_contains_oneL(self):
        r = _make_residue(oneL="W")
        assert "W" in repr(r)

    def test_repr_contains_rtype(self):
        r = _make_residue(rtype="Insertion")
        assert "Insertion" in repr(r)
