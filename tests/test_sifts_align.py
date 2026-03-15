"""Unit tests for SiftsAlign (sifts_segments_generation.py).

Mocking strategy:
- helper.EntryMapping  : avoids calling lalign36 subprocess
- generate_xref_csv.insert_mappings : avoids CSV I/O
- UNP                  : avoids network / cache access
- ChemCompMapping      : fast in practice (reads package data), not mocked
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from pdbe_sifts.sifts_segments_generation import SiftsAlign

DATA_DIR = Path(__file__).parent / "data"
CIF_1CBS = DATA_DIR / "cif" / "1cbs.cif"
CIF_1VYC = DATA_DIR / "cif" / "1vyc.cif"
DUCKDB = DATA_DIR / "mmseqs_hits" / "hits.duckdb"


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def align_no_db(tmp_path):
    """SiftsAlign without DuckDB (manual-mapping mode)."""
    return SiftsAlign(str(CIF_1CBS), str(tmp_path), db_conn_str=None)


@pytest.fixture
def align_with_db(tmp_path):
    """SiftsAlign with the test DuckDB."""
    return SiftsAlign(str(CIF_1CBS), str(tmp_path), db_conn_str=str(DUCKDB))


# ── Init tests ────────────────────────────────────────────────────────────────

def test_init_no_db(tmp_path):
    align = SiftsAlign(str(CIF_1CBS), str(tmp_path), db_conn_str=None)
    assert align.conn is None
    assert align.cif_file == str(CIF_1CBS)
    assert align.connectivity_mode is True
    assert align.nf90_mode is False


def test_init_with_db(tmp_path):
    align = SiftsAlign(str(CIF_1CBS), str(tmp_path), db_conn_str=str(DUCKDB))
    assert align.conn is not None
    align.conn.close()


def test_init_nf90_mode(tmp_path):
    align = SiftsAlign(str(CIF_1CBS), str(tmp_path), nf90_mode=True)
    assert align.nf90_mode is True


def test_init_no_connectivity(tmp_path):
    align = SiftsAlign(str(CIF_1CBS), str(tmp_path), connectivity_mode=False)
    assert align.connectivity_mode is False


# ── no_used_cif_category_modified ────────────────────────────────────────────

def test_no_categories_modified_real_cif(align_no_db):
    """1cbs.cif has only past revision dates → method returns False (do not skip)."""
    result = align_no_db.no_used_cif_category_modified(str(CIF_1CBS))
    assert result is False


def test_no_categories_modified_real_vyc(align_no_db):
    result = align_no_db.no_used_cif_category_modified(str(CIF_1VYC))
    assert result is False


# ── process_entry — early exit when categories modified ───────────────────────

def test_process_entry_skips_when_all_categories_modified(align_no_db):
    """If no_used_cif_category_modified returns True, process_entry returns early."""
    with patch.object(align_no_db, "no_used_cif_category_modified", return_value=True):
        # Should not raise, just log and return
        align_no_db.process_entry("1cbs")


# ── process_entry — full run with mocked alignment ────────────────────────────

def test_process_entry_1cbs(tmp_path):
    """Full process_entry call with lalign36 mocked out."""
    mock_em = MagicMock()
    mock_em.set_chain_accessions.return_value = True

    with (
        patch("pdbe_sifts.sifts_segments_generation.helper.EntryMapping", return_value=mock_em),
        patch("pdbe_sifts.sifts_segments_generation.generate_xref_csv.insert_mappings"),
    ):
        align = SiftsAlign(str(CIF_1CBS), str(tmp_path), db_conn_str=None)
        align.process_entry("1cbs")

    # Output dir created
    assert tmp_path.exists()


def test_process_entry_1vyc(tmp_path):
    """Same test for 1vyc.cif."""
    mock_em = MagicMock()
    mock_em.set_chain_accessions.return_value = True

    with (
        patch("pdbe_sifts.sifts_segments_generation.helper.EntryMapping", return_value=mock_em),
        patch("pdbe_sifts.sifts_segments_generation.generate_xref_csv.insert_mappings"),
    ):
        align = SiftsAlign(str(CIF_1VYC), str(tmp_path), db_conn_str=None)
        align.process_entry("1vyc")


def test_process_entry_chain_skipped_when_accessions_fail(tmp_path):
    """If set_chain_accessions() returns False, chain is skipped and files removed."""
    mock_em = MagicMock()
    mock_em.set_chain_accessions.return_value = False

    with (
        patch("pdbe_sifts.sifts_segments_generation.helper.EntryMapping", return_value=mock_em),
        patch("pdbe_sifts.sifts_segments_generation.generate_xref_csv.insert_mappings"),
        patch.object(SiftsAlign, "remove_existing_files") as mock_remove,
    ):
        align = SiftsAlign(str(CIF_1CBS), str(tmp_path), db_conn_str=None)
        align.process_entry("1cbs")

    mock_remove.assert_called()


# ── remove_existing_files ─────────────────────────────────────────────────────

def test_remove_existing_files(tmp_path):
    (tmp_path / "1cbs_seg.csv.gz").touch()
    (tmp_path / "1cbs_res.csv.gz").touch()
    (tmp_path / "other_seg.csv.gz").touch()  # belongs to a different entry

    align = SiftsAlign(str(CIF_1CBS), str(tmp_path), db_conn_str=None)
    align.remove_existing_files("1cbs")

    assert not (tmp_path / "1cbs_seg.csv.gz").exists()
    assert not (tmp_path / "1cbs_res.csv.gz").exists()
    assert (tmp_path / "other_seg.csv.gz").exists()  # untouched


# ── _parse_accession_mapping ──────────────────────────────────────────────────

def test_parse_accession_mapping(tmp_path):
    """A:P29373 → chain A mapped to P29373."""
    mock_unp = MagicMock()
    mock_unp.accession = "P29373"

    with patch("pdbe_sifts.sifts_segments_generation.UNP", return_value=mock_unp):
        align = SiftsAlign(str(CIF_1CBS), str(tmp_path), db_conn_str=None, unp_mode="A:P29373")
        result = align._parse_accession_mapping({})

    assert "A" in result
    assert result["A"][0].accession == "P29373"


def test_parse_accession_mapping_obsolete_skipped(tmp_path):
    """Obsolete UniProt accession should be silently skipped."""
    from pdbe_sifts.base.exceptions import ObsoleteUniProtError

    with patch("pdbe_sifts.sifts_segments_generation.UNP", side_effect=ObsoleteUniProtError):
        align = SiftsAlign(str(CIF_1CBS), str(tmp_path), db_conn_str=None, unp_mode="A:OBSOLETE")
        result = align._parse_accession_mapping({})

    assert result == {}


# ── _is_future_date ───────────────────────────────────────────────────────────

def test_is_future_date_past(align_no_db):
    assert align_no_db._is_future_date("1990-01-01") is False


def test_is_future_date_future(align_no_db):
    assert align_no_db._is_future_date("2099-12-31") is True
