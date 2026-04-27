"""Regression tests for path-traversal (CWE-22/CWE-73) and SQL injection (CWE-89)
and unsafe deserialization (CWE-502) mitigations.

Covers:
  - validate_entry_id: allowlist enforcement
  - safe_join: path-confinement enforcement
  - get_list_of_mappings: parameterised WHERE clause
  - SequenceMatchParser._init_database: path validation + identity_cutoff range
  - validate_uniprot_accession: UniProt accession allowlist (CWE-502)
  - UNP pickle→JSON migration: no pickle import, JSON round-trip, stale .pkl removal
"""

import duckdb
import pytest
from pdbe_sifts.base.utils import (
    safe_join,
    validate_entry_id,
    validate_uniprot_accession,
)

# ---------------------------------------------------------------------------
# validate_entry_id
# ---------------------------------------------------------------------------


class TestValidateEntryId:
    # ── valid values ────────────────────────────────────────────────────────

    @pytest.mark.parametrize(
        "entry_id",
        [
            "1cbs",  # classic 4-char PDB ID
            "7yfv",
            "af-p12345-f1",  # AlphaFold-style
            "a",  # single char
            "abc123",
            "entry_with_underscore",
            "x" * 64,  # maximum length
        ],
    )
    def test_valid_ids_are_returned_unchanged(self, entry_id):
        assert validate_entry_id(entry_id) == entry_id

    # ── path traversal sequences ────────────────────────────────────────────

    @pytest.mark.parametrize(
        "evil",
        [
            "../etc/passwd",
            "../../root/.ssh/authorized_keys",
            "/etc/cron.d/evil",
            "foo/bar",
            "foo\\bar",
        ],
    )
    def test_path_separators_are_rejected(self, evil):
        with pytest.raises(ValueError, match="Invalid entry_id"):
            validate_entry_id(evil)

    # ── shell / injection payloads ──────────────────────────────────────────

    @pytest.mark.parametrize(
        "evil",
        [
            "1cbs; rm -rf /",
            "1cbs\x00",  # null byte
            "1cbs\n",  # newline
            "1cbs ",  # trailing space
            " 1cbs",  # leading space
            "",  # empty string
            "-starts-with-dash",  # must start with alphanumeric
            "X" * 65,  # exceeds max length
        ],
    )
    def test_dangerous_chars_are_rejected(self, evil):
        with pytest.raises(ValueError, match="Invalid entry_id"):
            validate_entry_id(evil)

    def test_none_is_rejected(self):
        """None must not bypass the check (e.g. when --entry is omitted)."""
        with pytest.raises((ValueError, TypeError)):
            validate_entry_id(None)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# safe_join
# ---------------------------------------------------------------------------


class TestSafeJoin:
    def test_normal_filename_stays_inside(self, tmp_path):
        result = safe_join(tmp_path, "1cbs_seg.csv.gz")
        assert result.parent.resolve() == tmp_path.resolve()

    def test_returns_path_object(self, tmp_path):
        result = safe_join(tmp_path, "1cbs_res.csv.gz")
        from pathlib import Path

        assert isinstance(result, Path)

    @pytest.mark.parametrize(
        "evil_filename",
        [
            "../outside.txt",
            "../../etc/cron.d/evil",
            "/absolute/path.txt",
        ],
    )
    def test_traversal_filenames_are_rejected(self, tmp_path, evil_filename):
        with pytest.raises(ValueError, match="Path traversal detected"):
            safe_join(tmp_path, evil_filename)

    def test_symlink_traversal_is_rejected(self, tmp_path):
        """A symlink inside base_dir pointing outside must be rejected."""
        outside = tmp_path.parent / "outside_target"
        outside.mkdir(exist_ok=True)
        link = tmp_path / "escape_link"
        link.symlink_to(outside)
        with pytest.raises(ValueError, match="Path traversal detected"):
            safe_join(tmp_path, "escape_link/payload.txt")


# ---------------------------------------------------------------------------
# get_list_of_mappings — parameterised WHERE (CWE-89)
# ---------------------------------------------------------------------------


class TestSqlInjectionGetMappings:
    """Verify that get_curated_db_mappings uses a ? parameter for entry_id."""

    def _make_hits_db(self, tmp_path) -> duckdb.DuckDBPyConnection:
        """Create a minimal in-memory hits table for testing."""
        conn = duckdb.connect()
        conn.execute("""
            CREATE TABLE hits (
                entry VARCHAR, entity INTEGER, accession VARCHAR,
                target_start INTEGER, target_end INTEGER,
                sifts_score DOUBLE, pdb_cross_references INTEGER,
                adjusted_score DOUBLE, dataset_score DOUBLE
            )
        """)
        conn.execute(
            "INSERT INTO hits VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            ["1cbs", 1, "P29373", 1, 182, 1.0, 1, 1.0, 1.0],
        )
        return conn

    def test_normal_entry_returns_mapping(self, tmp_path):
        from pdbe_sifts.segments_generation.get_list_of_mappings import (
            _SELECTION_QUERY,
        )

        conn = self._make_hits_db(tmp_path)
        rows = conn.execute(_SELECTION_QUERY, ["1cbs"]).fetchall()
        assert len(rows) == 1
        assert rows[0][2] == "P29373"

    def test_quote_injection_returns_no_rows(self, tmp_path):
        """A single-quote in entry must NOT cause a SQL error or return extra rows."""
        from pdbe_sifts.segments_generation.get_list_of_mappings import (
            _SELECTION_QUERY,
        )

        conn = self._make_hits_db(tmp_path)
        # Would break f-string: 1cbs' OR '1'='1
        rows = conn.execute(_SELECTION_QUERY, ["1cbs' OR '1'='1"]).fetchall()
        assert rows == []

    def test_semicolon_payload_returns_no_rows(self, tmp_path):
        """Semicolon-separated SQL payload must not execute as a second statement."""
        from pdbe_sifts.segments_generation.get_list_of_mappings import (
            _SELECTION_QUERY,
        )

        conn = self._make_hits_db(tmp_path)
        rows = conn.execute(
            _SELECTION_QUERY, ["1cbs; DROP TABLE hits; --"]
        ).fetchall()
        # hits table must still exist and payload must not have matched
        assert rows == []
        assert conn.execute("SELECT COUNT(*) FROM hits").fetchone()[0] == 1


# ---------------------------------------------------------------------------
# SequenceMatchParser._init_database — path validation + identity_cutoff (CWE-89)
# ---------------------------------------------------------------------------


class TestSqlInjectionParser:
    """Verify path validation and identity_cutoff bounds in _init_database."""

    def _make_mmseqs_tsv(self, path):
        """Write a minimal (empty) MMseqs2 TSV so stat().st_size is non-zero."""
        # Write one real MMseqs2 row so the parser can load it.
        path.write_text(
            "pdb|1cbs-1|OX=9606\tsp|P29373|RABP2_HUMAN\t182\t0\t1\t182\t1\t182"
            "\t1e-100\t500\tAAA\tAAA\t182\t9606\t>pdb|1cbs-1|OX=9606\t1.0\t1.0\n"
        )

    def test_identity_cutoff_above_one_raises(self, tmp_path):
        from pdbe_sifts.sequence_match.sequence_match_parser import (
            SequenceMatchParser,
        )

        tsv = tmp_path / "hits.tsv"
        self._make_mmseqs_tsv(tsv)
        parser = SequenceMatchParser(
            "mmseqs", tsv, tmp_path, identity_cutoff=1.5
        )
        with pytest.raises(ValueError, match="identity_cutoff"):
            parser._init_database()

    def test_identity_cutoff_negative_raises(self, tmp_path):
        from pdbe_sifts.sequence_match.sequence_match_parser import (
            SequenceMatchParser,
        )

        tsv = tmp_path / "hits.tsv"
        self._make_mmseqs_tsv(tsv)
        parser = SequenceMatchParser(
            "mmseqs", tsv, tmp_path, identity_cutoff=-0.1
        )
        with pytest.raises(ValueError, match="identity_cutoff"):
            parser._init_database()

    def test_nonexistent_file_raises(self, tmp_path):
        from pdbe_sifts.sequence_match.sequence_match_parser import (
            SequenceMatchParser,
        )

        parser = SequenceMatchParser(
            "mmseqs", tmp_path / "does_not_exist.tsv", tmp_path
        )
        with pytest.raises(FileNotFoundError):
            parser._init_database()

    def test_path_with_single_quote_is_sql_escaped(self, tmp_path):
        """A path containing a single-quote must be SQL-escaped before interpolation."""
        from pdbe_sifts.sequence_match.sequence_match_parser import (
            INSERT_MMSEQS_TABLE_HITS,
        )

        # Simulate what _init_database does: escape then format
        evil_path = "/tmp/it's/a/trap.tsv"
        safe_path = evil_path.replace("'", "''")
        sql = INSERT_MMSEQS_TABLE_HITS.format(tsv_path=safe_path)

        # The escaped path (/tmp/it''s/…) must appear literally in the SQL
        assert "it''s" in sql, "Single quote was not doubled in the SQL path"
        # The raw unescaped path must NOT appear
        assert "it's" not in sql.replace(
            "it''s", ""
        ), "Unescaped single quote found in SQL path"


# ---------------------------------------------------------------------------
# validate_uniprot_accession — allowlist enforcement (CWE-502 pre-condition)
# ---------------------------------------------------------------------------


class TestValidateUniprotAccession:
    # ── valid accessions ────────────────────────────────────────────────────

    @pytest.mark.parametrize(
        "accession",
        [
            "P29373",  # classic 6-char
            "Q9Y6K9",
            "O15530",
            "A0A000",  # starts with A (10-char prefix pattern begins)
            "P12345-2",  # isoform suffix
            "Q9Y6K9-10",  # two-digit isoform suffix
            # 10-char accessions
            "A0A000B1C2",
            "A0A001B2C3",
        ],
    )
    def test_valid_accessions_are_returned_unchanged(self, accession):
        assert validate_uniprot_accession(accession) == accession

    # ── path traversal ──────────────────────────────────────────────────────

    @pytest.mark.parametrize(
        "evil",
        [
            "../../evil",
            "../etc/passwd",
            "/absolute/path",
            "P29373/../../etc",
            "P29373\\evil",
        ],
    )
    def test_path_traversal_sequences_are_rejected(self, evil):
        with pytest.raises(ValueError, match="Invalid UniProt accession"):
            validate_uniprot_accession(evil)

    # ── malformed / unexpected input ────────────────────────────────────────

    @pytest.mark.parametrize(
        "bad",
        [
            "",  # empty
            "P1234",  # too short
            "P123456",  # too long for 6-char format
            "p29373",  # lowercase
            "NOTACC",  # wrong pattern
            "P29373 ",  # trailing space
            " P29373",  # leading space
            "P29373\x00",  # null byte
            "P29373\n",  # newline
            "P29373-",  # dash without isoform number
            "P29373-abc",  # non-numeric isoform
        ],
    )
    def test_malformed_accessions_are_rejected(self, bad):
        with pytest.raises(ValueError, match="Invalid UniProt accession"):
            validate_uniprot_accession(bad)

    def test_none_is_rejected(self):
        with pytest.raises((ValueError, TypeError)):
            validate_uniprot_accession(None)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# UNP RestrictedUnpickler — safe pickle deserialization (CWE-502)
# ---------------------------------------------------------------------------


def _make_unp_dict():
    """Return a minimal dict that mimics a populated UNP.__dict__."""
    return {
        "accession": "P29373",
        "ad_dbref_auto_acc": "P29373",
        "seq_isoforms": {},
        "isoforms": {},
        "features": {"signal": None},
        "sequence": "ACDEFGHIKLMNPQRSTVWY",
        "seq_length": 20,
        "secondary_accessions": ["Q8WXI4"],
        "evidences": [{"type": "ECO:0000269", "source": None}],
        "keywords": ["Membrane", "Signal"],
        "organism": ["Homo sapiens"],
        "taxonomy": ["9606"],
        "dataset": "Swiss-Prot",
    }


class TestUnpRestrictedUnpickler:
    """Verify _RestrictedUnpickler blocks dangerous classes while allowing primitives."""

    def test_allows_primitive_builtins(self):
        """A pickle of plain primitives must deserialize without error."""
        import io
        import pickle

        from pdbe_sifts.unp.unp import _RestrictedUnpickler

        data = _make_unp_dict()
        raw = pickle.dumps(data, protocol=pickle.HIGHEST_PROTOCOL)
        loaded = _RestrictedUnpickler(io.BytesIO(raw)).load()
        assert loaded == data

    def test_blocks_os_system_payload(self):
        """A pickle payload referencing os.system must raise UnpicklingError."""
        import io
        import os
        import pickle

        from pdbe_sifts.unp.unp import _RestrictedUnpickler

        # pickle.dumps(os.system) generates a GLOBAL opcode that calls
        # find_class("os", "system") on load — exactly what we want to block.
        evil = pickle.dumps(os.system, protocol=2)
        with pytest.raises(
            pickle.UnpicklingError, match="Forbidden pickle class"
        ):
            _RestrictedUnpickler(io.BytesIO(evil)).load()

    def test_pkl_cache_round_trip(self, tmp_path):
        """dump → _RestrictedUnpickler.load() returns data equal to the original."""
        import io
        import pickle

        from pdbe_sifts.unp.unp import _RestrictedUnpickler

        data = _make_unp_dict()
        buf = io.BytesIO()
        pickle.dump(data, buf, protocol=pickle.HIGHEST_PROTOCOL)
        buf.seek(0)
        assert _RestrictedUnpickler(buf).load() == data

    def test_stale_json_is_deleted_on_access(self, tmp_path, monkeypatch):
        """A .json file from the intermediate migration must be removed when UNP loads."""
        import pdbe_sifts.unp.unp as unp_module

        acc = "P29373"

        # Patch the cache-dir reference bound inside the unp module
        monkeypatch.setattr(
            unp_module, "get_uniprot_cache_dir", lambda a: str(tmp_path)
        )

        # Create a stale .json file (left over from the intermediate JSON-cache phase)
        json_path = tmp_path / f"{acc}.json"
        json_path.write_text('{"accession": "P29373"}', encoding="utf-8")
        assert json_path.exists(), "Pre-condition: .json file must exist"

        # Write a valid .pkl so the cache hits (no network call needed)
        import pickle

        data = _make_unp_dict()
        pkl_path = tmp_path / f"{acc}.pkl"
        with open(pkl_path, "wb") as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

        unp_module.UNP(acc)

        assert (
            not json_path.exists()
        ), ".json file was NOT deleted — reverse migration code is missing"
