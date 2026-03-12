import pytest
from pathlib import Path
from unittest.mock import MagicMock

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def cif_list_file(tmp_path):
    p = tmp_path / "cif_path.txt"
    cif_dir = DATA_DIR / "cif"
    p.write_text(f"{cif_dir / '1cbs.cif'}\n{cif_dir / '1vyc.cif'}\n")
    return p


@pytest.fixture
def mock_ncbi():
    lineages_leaf_to_root = {
        10665:   [10665, 10663, 1198136, 2946170, 3420545, 2731619, 2731618, 2731360, 2731341, 10239, 1],
        697290:  [697290, 10665, 10663, 1198136, 2946170, 3420545, 2731619, 2731618, 2731360, 2731341, 10239, 1],
        2750851: [2750851, 2844216, 10663, 1198136, 2946170, 3420545, 2731619, 2731618, 2731360, 2731341, 10239, 1],
        2696339: [2696339, 2844193, 3430963, 1913652, 1198136, 2946170, 3420545, 2731619, 2731618, 2731360, 2731341, 10239, 1],
        3044333: [3044333, 12333, 10239, 1],
    }
    mock = MagicMock()
    mock.get_lineage.side_effect = lambda taxid: list(reversed(lineages_leaf_to_root[taxid]))
    return mock


@pytest.fixture
def unp_csv(tmp_path):
    p = tmp_path / "unp.csv"
    p.write_text(
        ",accession,provenance,pdb_xref,annotation_score\n"
        "0,P29373,Swiss-Prot,81,5.0\n"
        "1,P83346,Swiss-Prot,2,3.0\n"
    )
    return p
