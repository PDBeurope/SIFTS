from pathlib import Path

import pytest
from pdbe_sifts.mmcif import extract_entities

DATA_DIR = Path(__file__).parent / "data" / "cif"


def test_returns_sequence_and_taxid(tmp_path):
    result = extract_entities(DATA_DIR / "1cbs.cif")
    assert "1" in result
    sequence, tax_id = result["1"]
    assert (
        sequence
        == "PNFSGNWKIIRSENFEELLKVLGVNVMLRKIAVAAASKPAVEIKQEGDTFYIKTSTTVRTTEINFKVGEEFEEQTVDGRPCKSLVKWESENKMVCEQKLLKGEGPKTSWTRELTNDGELILTMTADDVVCTRVYVRE"
    )
    assert tax_id == 9606


def test_sequence_length(tmp_path):
    result = extract_entities(DATA_DIR / "1cbs.cif")
    sequence, _ = result["1"]
    assert len(sequence) == 137


def test_only_polypeptide_entities_returned():
    result = extract_entities(DATA_DIR / "1cbs.cif")
    assert all(
        isinstance(seq, str) and len(seq) > 0 for seq, _ in result.values()
    )


def test_second_entry(tmp_path):
    result = extract_entities(DATA_DIR / "1vyc.cif")
    assert "1" in result
    sequence, tax_id = result["1"]
    assert (
        sequence
        == "RKCLIKYSQANESSKTCPSGQLLCLKKWEIGNPSGKEVKRGCVATCPKPWKNEIIQCCAKDKCNA"
    )
    assert tax_id == 92438


def test_file_not_found():
    with pytest.raises(FileNotFoundError):
        extract_entities("nonexistent.cif")


def test_returns_dict():
    result = extract_entities(DATA_DIR / "1cbs.cif")
    assert isinstance(result, dict)
    assert len(result) > 0
