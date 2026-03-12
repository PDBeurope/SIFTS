import duckdb
import pytest


def _make_hits(rows: list[tuple]) -> duckdb.DuckDBPyConnection:
    conn = duckdb.connect()
    conn.execute("""
        CREATE TABLE hits (
            identity      DOUBLE,
            coverage      DOUBLE,
            mismatch      INTEGER,
            query_len     INTEGER,
            adjusted_score DOUBLE,
            tax_score      DOUBLE,
            dataset_score  DOUBLE,
            sifts_score    DOUBLE
        )
    """)
    conn.executemany("INSERT INTO hits VALUES (?, ?, ?, ?, ?, ?, ?, ?)", rows)
    return conn


def _make_accession_info(rows: list[tuple]) -> duckdb.DuckDBPyConnection:
    conn = duckdb.connect()
    conn.execute("""
        CREATE TABLE accession_info (
            provenance       VARCHAR,
            annotation_score DOUBLE,
            dataset_score    DOUBLE
        )
    """)
    conn.executemany("INSERT INTO accession_info VALUES (?, ?, NULL)", rows)
    return conn


_ADJUSTED_SQL = """
    UPDATE hits
    SET adjusted_score =
        CASE
            WHEN identity IS NULL
            OR coverage IS NULL
            OR mismatch IS NULL
            OR query_len IS NULL
            OR query_len = 0
            THEN 0
            ELSE (identity * coverage) * (1 - (mismatch::DOUBLE / query_len)) * 1000
        END
"""

_DATASET_SQL = """
    UPDATE accession_info
    SET dataset_score =
        CASE
            WHEN annotation_score IS NULL THEN 0
            WHEN provenance = 'Swiss-Prot'
                THEN 20 * annotation_score
            ELSE
                -50 + 10 * annotation_score
        END
"""

_SIFTS_SQL = """
    UPDATE hits
    SET sifts_score =
        COALESCE(adjusted_score, 0)
        + COALESCE(tax_score, 0)
        + COALESCE(dataset_score, 0)
"""


# ---------------------------------------------------------------------------
# adjusted_score
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("identity, coverage, mismatch, query_len, expected", [
    (1.0,  1.0,  0,   137, 1000.0),
    (1.0,  1.0,  0,    65, 1000.0),
    (0.95, 0.90, 5,   100,  0.95 * 0.90 * (1 - 5 / 100) * 1000),
    (0.50, 0.50, 10,   50,  0.50 * 0.50 * (1 - 10 / 50) * 1000),
    (1.0,  1.0,  0,     0,  0.0),
    (None, 1.0,  0,   100,  0.0),
    (1.0,  None, 0,   100,  0.0),
    (1.0,  1.0,  None, 100, 0.0),
])
def test_adjusted_score_formula(identity, coverage, mismatch, query_len, expected):
    conn = _make_hits([(identity, coverage, mismatch, query_len, None, None, None, None)])
    conn.execute(_ADJUSTED_SQL)
    result = conn.execute("SELECT adjusted_score FROM hits").fetchone()[0]
    assert result == pytest.approx(expected)


# ---------------------------------------------------------------------------
# dataset_score
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("provenance, annotation_score, expected", [
    ("Swiss-Prot", 5.0,  100.0),
    ("Swiss-Prot", 3.0,   60.0),
    ("Swiss-Prot", 1.0,   20.0),
    ("TrEMBL",     5.0,    0.0),   # -50 + 10*5
    ("TrEMBL",     6.0,   10.0),   # -50 + 10*6
    ("TrEMBL",     1.0,  -40.0),   # -50 + 10*1
    ("Swiss-Prot", None,   0.0),
    (None,         5.0,  -50 + 10 * 5),   # provenance None → ELSE branch
])
def test_dataset_score_formula(provenance, annotation_score, expected):
    conn = _make_accession_info([(provenance, annotation_score)])
    conn.execute(_DATASET_SQL)
    result = conn.execute("SELECT dataset_score FROM accession_info").fetchone()[0]
    assert result == pytest.approx(expected)


# ---------------------------------------------------------------------------
# sifts_score = adjusted + tax + dataset  (COALESCE nulls to 0)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("adjusted, tax, dataset, expected", [
    (1000.0, 200.0, 100.0, 1300.0),
    (1000.0, 200.0,  60.0, 1260.0),
    (1000.0, None,  100.0, 1100.0),   # tax NULL → 0
    (None,   200.0, 100.0,  300.0),   # adjusted NULL → 0
    (None,   None,  None,     0.0),   # all NULL → 0
    (500.0,  50.0,  -40.0,  510.0),   # dataset négatif (TrEMBL bas score)
])
def test_sifts_score_formula(adjusted, tax, dataset, expected):
    conn = duckdb.connect()
    conn.execute("""
        CREATE TABLE hits (
            adjusted_score DOUBLE,
            tax_score      DOUBLE,
            dataset_score  DOUBLE,
            sifts_score    DOUBLE
        )
    """)
    conn.execute("INSERT INTO hits VALUES (?, ?, ?, NULL)", [adjusted, tax, dataset])
    conn.execute(_SIFTS_SQL)
    result = conn.execute("SELECT sifts_score FROM hits").fetchone()[0]
    assert result == pytest.approx(expected)
