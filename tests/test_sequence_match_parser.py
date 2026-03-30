from pathlib import Path

import duckdb
import pytest
from pdbe_sifts.sequence_match.sequence_match_parser import SequenceMatchParser

DATA_DIR = Path(__file__).parent / "data"
BLASTP_TSV = DATA_DIR / "blastp_hits" / "hits_tmp.tsv"
MMSEQS_TSV = DATA_DIR / "mmseqs_hits" / "hits_tmp.tsv"


def make_parser(fmt, tsv, tmp_path, unp_csv=None):
    return SequenceMatchParser(
        format=fmt,
        result_file_path=tsv,
        out_dir=tmp_path,
        unp_csv=unp_csv,
        max_workers=1,
        identity_cutoff=0.9,
    )


def query(db_path, sql, params=None):
    conn = duckdb.connect(str(db_path), read_only=True)
    result = conn.execute(sql, params or []).fetchall()
    conn.close()
    return result


# ---------------------------------------------------------------------------
# _init_database
# ---------------------------------------------------------------------------


def test_blastp_loads_two_rows(tmp_path):
    p = make_parser("blastp", BLASTP_TSV, tmp_path)
    p._init_database()
    count = query(p.db_path, "SELECT COUNT(*) FROM hits")[0][0]
    assert count == 2


def test_mmseqs_loads_two_rows(tmp_path):
    p = make_parser("mmseqs", MMSEQS_TSV, tmp_path)
    p._init_database()
    count = query(p.db_path, "SELECT COUNT(*) FROM hits")[0][0]
    assert count == 2


def test_blastp_entry_entity_parsed(tmp_path):
    p = make_parser("blastp", BLASTP_TSV, tmp_path)
    p._init_database()
    rows = query(
        p.db_path, "SELECT entry, entity, accession FROM hits ORDER BY entry"
    )
    assert ("1cbs", "1", "P29373") in rows
    assert ("1vyc", "1", "P83346") in rows


def test_mmseqs_entry_entity_parsed(tmp_path):
    p = make_parser("mmseqs", MMSEQS_TSV, tmp_path)
    p._init_database()
    rows = query(
        p.db_path, "SELECT entry, entity, accession FROM hits ORDER BY entry"
    )
    assert ("1cbs", "1", "P29373") in rows
    assert ("1vyc", "1", "P83346") in rows


def test_blastp_identity_coverage(tmp_path):
    p = make_parser("blastp", BLASTP_TSV, tmp_path)
    p._init_database()
    rows = query(p.db_path, "SELECT identity, coverage FROM hits")
    for identity, coverage in rows:
        assert identity == pytest.approx(1.0)
        assert coverage == pytest.approx(1.0)


def test_blastp_tax_ids_parsed(tmp_path):
    p = make_parser("blastp", BLASTP_TSV, tmp_path)
    p._init_database()
    rows = {
        r[0]: (r[1], r[2])
        for r in query(
            p.db_path, "SELECT entry, query_tax_id, target_tax_id FROM hits"
        )
    }
    assert rows["1cbs"] == (9606, 9606)
    assert rows["1vyc"] == (92438, 92438)


# ---------------------------------------------------------------------------
# compute_adjusted_score
# ---------------------------------------------------------------------------


def test_adjusted_score_blastp(tmp_path):
    p = make_parser("blastp", BLASTP_TSV, tmp_path)
    p._init_database()
    p.compute_adjusted_score()
    rows = query(p.db_path, "SELECT adjusted_score FROM hits")
    for (score,) in rows:
        assert score == pytest.approx(1000.0)


def test_adjusted_score_mmseqs(tmp_path):
    p = make_parser("mmseqs", MMSEQS_TSV, tmp_path)
    p._init_database()
    p.compute_adjusted_score()
    rows = query(p.db_path, "SELECT adjusted_score FROM hits")
    for (score,) in rows:
        assert score == pytest.approx(1000.0)


# ---------------------------------------------------------------------------
# compute_dataset_score (requires unp_csv fixture)
# ---------------------------------------------------------------------------


def test_dataset_score_with_unp_csv(tmp_path, unp_csv):
    p = make_parser("blastp", BLASTP_TSV, tmp_path, unp_csv=unp_csv)
    p._init_database()
    p.compute_dataset_score()
    rows = {
        r[0]: r[1]
        for r in query(p.db_path, "SELECT accession, dataset_score FROM hits")
    }
    assert rows["P29373"] == pytest.approx(100.0)
    assert rows["P83346"] == pytest.approx(60.0)


# ---------------------------------------------------------------------------
# compute_tax_score
# Same taxid in query and target → get_tax_weight returns 200 immediately
# without touching the ete4 NCBI database.
# ---------------------------------------------------------------------------


def test_tax_score_same_taxid(tmp_path):
    p = make_parser("blastp", BLASTP_TSV, tmp_path)
    p._init_database()
    p.compute_tax_score()
    rows = query(p.db_path, "SELECT tax_score FROM hits")
    for (score,) in rows:
        assert score == pytest.approx(200.0)


# ---------------------------------------------------------------------------
# compute_sifts_score (full pipeline without search)
# adjusted=1000, tax=200, dataset: P29373→100, P83346→60
# ---------------------------------------------------------------------------


def test_sifts_score_full(tmp_path, unp_csv):
    p = make_parser("blastp", BLASTP_TSV, tmp_path, unp_csv=unp_csv)
    p._init_database()
    p.compute_adjusted_score()
    p.compute_tax_score()
    p.compute_dataset_score()
    p.compute_sifts_score()

    rows = {
        r[0]: r[1]
        for r in query(p.db_path, "SELECT accession, sifts_score FROM hits")
    }
    assert rows["P29373"] == pytest.approx(1300.0)
    assert rows["P83346"] == pytest.approx(1260.0)


# ---------------------------------------------------------------------------
# compute_rank
# ---------------------------------------------------------------------------


def test_rank_assigned(tmp_path, unp_csv):
    p = make_parser("blastp", BLASTP_TSV, tmp_path, unp_csv=unp_csv)
    p._init_database()
    p.compute_adjusted_score()
    p.compute_tax_score()
    p.compute_dataset_score()
    p.compute_sifts_score()
    p.compute_rank()

    rows = query(p.db_path, "SELECT hit_rank FROM hits")
    assert all(r[0] == 1 for r in rows)
