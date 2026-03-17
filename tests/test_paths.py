"""Tests for pdbe_sifts.base.paths — using base_dir to bypass conf."""

from pdbe_sifts.base.paths import ccd_cache_path, uniprot_cache_dir


class TestUniprotCacheDir:
    def test_structure(self):
        path = uniprot_cache_dir("P00520", base_dir="/cache/uniprot")
        assert path == "/cache/uniprot/P/P00520"

    def test_first_char_subdir(self):
        path = uniprot_cache_dir("Q9Y6K9", base_dir="/tmp")
        assert "/Q/" in path

    def test_accession_in_path(self):
        path = uniprot_cache_dir("A0A000", base_dir="/tmp")
        assert "A0A000" in path


class TestCcdCachePath:
    def test_structure(self):
        path = ccd_cache_path("ALA", base_dir="/cache/ccd")
        assert path == "/cache/ccd/A/ALA/ALA.cif"

    def test_lowercase_is_uppercased(self):
        path = ccd_cache_path("ala", base_dir="/tmp")
        assert "ALA" in path
        assert "ala" not in path

    def test_ends_with_cif(self):
        path = ccd_cache_path("GLY", base_dir="/tmp")
        assert path.endswith(".cif")

    def test_first_char_subdir(self):
        path = ccd_cache_path("MET", base_dir="/tmp")
        assert "/M/" in path
