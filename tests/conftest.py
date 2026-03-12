import pytest
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def cif_list_file(tmp_path):
    p = tmp_path / "cif_path.txt"
    cif_dir = DATA_DIR / "cif"
    p.write_text(f"{cif_dir / '1cbs.cif'}\n{cif_dir / '1vyc.cif'}\n")
    return p


@pytest.fixture
def unp_csv(tmp_path):
    p = tmp_path / "unp.csv"
    p.write_text(
        ",accession,provenance,pdb_xref,annotation_score\n"
        "0,P29373,Swiss-Prot,81,5.0\n"
        "1,P83346,Swiss-Prot,2,3.0\n"
    )
    return p
