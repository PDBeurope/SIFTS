"""Tests for pdbe_sifts.segments_generation.connectivity.ccd_parser."""

from pathlib import Path
from unittest.mock import patch

import pytest
from pdbe_sifts.segments_generation.connectivity.ccd_parser import (
    CcdFile,
    ccd_dict_builder,
    get_ccd_file,
)

CCD_FIXTURE = Path(__file__).parent / "data" / "ccd" / "ALA.cif"


class TestCcdDictBuilder:
    def test_returns_n_and_c_keys(self):
        ccd_ters = {
            "n_ter": (0, {"atom_id": "N", "pdbx_ordinal": "1"}),
            "c_ter": (4, {"atom_id": "C", "pdbx_ordinal": "5"}),
        }
        result = ccd_dict_builder(ccd_ters)
        assert "N" in result
        assert "C" in result

    def test_n_ter_atom_id(self):
        ccd_ters = {
            "n_ter": (0, {"atom_id": "N", "pdbx_ordinal": "1"}),
            "c_ter": (4, {"atom_id": "C", "pdbx_ordinal": "5"}),
        }
        result = ccd_dict_builder(ccd_ters)
        assert result["N"] == ("N", "1")

    def test_c_ter_atom_id(self):
        ccd_ters = {
            "n_ter": (0, {"atom_id": "N", "pdbx_ordinal": "1"}),
            "c_ter": (4, {"atom_id": "C", "pdbx_ordinal": "5"}),
        }
        result = ccd_dict_builder(ccd_ters)
        assert result["C"] == ("C", "5")

    def test_missing_key_returns_zero(self):
        # Missing c_ter → KeyError → returns 0
        result = ccd_dict_builder({"n_ter": (0, {"atom_id": "N"})})
        assert result == 0

    def test_empty_dict_returns_zero(self):
        result = ccd_dict_builder({})
        assert result == 0


class TestGetCcdFile:
    def test_returns_existing_path_without_download(self, tmp_path):
        # Pre-create a file in the expected cache location
        code = "ALA"
        ccd_dir = tmp_path / "A" / "ALA"
        ccd_dir.mkdir(parents=True)
        ccd_file = ccd_dir / "ALA.cif"
        ccd_file.write_text("data_ALA")

        with patch(
            "pdbe_sifts.segments_generation.connectivity.ccd_parser.ccd_cache_path",
            return_value=str(ccd_file),
        ):
            result = get_ccd_file(code)
            assert result == str(ccd_file)

    def test_triggers_download_when_missing(self, tmp_path):
        code = "GLY"
        ccd_file = tmp_path / "G" / "GLY" / "GLY.cif"

        with (
            patch(
                "pdbe_sifts.segments_generation.connectivity.ccd_parser.ccd_cache_path",
                return_value=str(ccd_file),
            ),
            patch(
                "pdbe_sifts.segments_generation.connectivity.ccd_parser.download_file_from_url"
            ) as mock_dl,
        ):
            # Simulate download writing the file
            def write_file(url, path):
                Path(path).parent.mkdir(parents=True, exist_ok=True)
                Path(path).write_text("data_GLY")

            mock_dl.side_effect = write_file
            result = get_ccd_file(code)
            assert mock_dl.called
            assert result == str(ccd_file)


@pytest.mark.skipif(
    not CCD_FIXTURE.exists(), reason="ALA.cif fixture not available"
)
class TestCcdFileWithFixture:
    @pytest.fixture
    def ala_ccd(self, tmp_path):
        ccd_path = tmp_path / "A" / "ALA" / "ALA.cif"
        ccd_path.parent.mkdir(parents=True)
        import shutil

        shutil.copy(CCD_FIXTURE, ccd_path)
        with patch(
            "pdbe_sifts.segments_generation.connectivity.ccd_parser.ccd_cache_path",
            return_value=str(ccd_path),
        ):
            yield CcdFile("ALA")

    def test_parse_ccd_populates_atoms(self, ala_ccd):
        ala_ccd.parse_ccd()
        assert ala_ccd.atoms is not None
        assert len(ala_ccd.atoms) > 0

    def test_extract_ters_returns_n_and_c(self, ala_ccd):
        ala_ccd.parse_ccd()
        ters = ala_ccd.extract_ters()
        assert "n_ter" in ters
        assert "c_ter" in ters

    def test_n_ter_is_nitrogen(self, ala_ccd):
        ala_ccd.parse_ccd()
        ters = ala_ccd.extract_ters()
        _, atom = ters["n_ter"]
        assert atom["type_symbol"] == "N"

    def test_c_ter_is_carbon(self, ala_ccd):
        ala_ccd.parse_ccd()
        ters = ala_ccd.extract_ters()
        _, atom = ters["c_ter"]
        assert atom["type_symbol"] == "C"

    def test_process_returns_dict(self, ala_ccd):
        result = ala_ccd.process()
        assert isinstance(result, dict)
        assert "N" in result
        assert "C" in result
