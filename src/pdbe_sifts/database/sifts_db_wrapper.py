from pathlib import Path

from pdbe_sifts.base.log import logger


def segment_to_dict(seg):
    """Serialise a SIFTS segment object to a plain dictionary.

    Args:
        seg: A segment object exposing the SIFTS segment attributes
            (``entry_id``, ``accession``, ``identity``, etc.).

    Returns:
        A ``dict`` mapping column names to their scalar values, ready for
        insertion into a pandas DataFrame or DuckDB table.
    """
    return {
        "entry_id": seg.entry_id,
        "entity_id": int(seg.entity_id),
        "id": seg.id,
        "auth_asym_id": seg.auth_asym_id,
        "struct_asym_id": seg.struct_asym_id,
        "accession": seg.accession,
        "name": seg.name,
        "seq_version": seg.seq_version,
        "unp_start": seg.unp_start,
        "pdb_start": seg.pdb_start,
        "unp_end": seg.unp_end,
        "pdb_end": seg.pdb_end,
        "auth_start": seg.auth_start,
        "auth_start_icode": seg.auth_start_icode,
        "auth_end": seg.auth_end,
        "auth_end_icode": seg.auth_end_icode,
        "conflicts": seg.conflicts,
        "modifications": seg.modifications,
        "unp_alignment": seg.unp_alignment,
        "pdb_alignment": seg.pdb_alignment,
        "identity": seg.identity,
        "score": seg.score,
        "best_mapping": bool(seg.best_mapping),
        "canonical_acc": bool(seg.canonical_acc),
        "reference_acc": seg.reference_acc,
        "chimera": bool(seg.chimera),
    }


def residue_to_dict(r):
    """Serialise a SIFTS residue object to a plain dictionary.

    Args:
        r: A residue object exposing the SIFTS residue attributes
            (``entry_id``, ``auth_seq_id``, ``unp_seq_id``, etc.).

    Returns:
        A ``dict`` mapping column names to their scalar values, ready for
        insertion into a pandas DataFrame or DuckDB table.
    """
    return {
        "entry_id": r.entry_id,
        "entity_id": int(r.entity_id),
        "id": r.id,
        "auth_asym_id": r.auth_asym_id,
        "struct_asym_id": r.struct_asym_id,
        "unp_segment_id": r.unp_segment_id,
        "auth_seq_id": r.auth_seq_id,
        "auth_seq_id_ins_code": r.auth_seq_id_ins_code,
        "pdb_seq_id": r.pdb_seq_id,
        "unp_seq_id": r.unp_seq_id,
        "observed": r.observed,
        "dbentry_id": r.dbentry_id,
        "accession": r.accession,
        "name": r.name,
        "type": r.type,
        "unp_one_letter_code": r.unp_one_letter_code,
        "pdb_one_letter_code": r.pdb_one_letter_code,
        "chem_comp_id": r.chem_comp_id,
        "mh_id": r.mh_id,
        "tax_id": r.tax_id,
        "canonical_acc": bool(r.canonical_acc),
        "reference_acc": r.reference_acc,
        "best_mapping": (
            bool(r.best_mapping) if r.best_mapping is not None else None
        ),
        "residue_id": r.residue_id,
    }


class SiftsDB:
    """Thin wrapper around a DuckDB connection for SIFTS cross-reference tables.

    On construction the two core tables (``sifts_xref_residue`` and
    ``sifts_xref_segment``) are created if they do not already exist.
    """

    def __init__(self, conn):
        """Initialise the wrapper and ensure both SIFTS tables exist.

        Args:
            conn: An open :class:`duckdb.DuckDBPyConnection` instance.
        """
        self.conn = conn
        self.create_sifts_xref_residue_table()
        self.create_sifts_xref_segment_table()

    def create_sifts_xref_residue_table(self):
        """Create the ``sifts_xref_residue`` table if it does not already exist."""
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS sifts_xref_residue (
                entry_id                VARCHAR NOT NULL,
                entity_id               INTEGER NOT NULL,
                id                      INTEGER NOT NULL,
                auth_asym_id            VARCHAR,
                struct_asym_id          VARCHAR NOT NULL,
                unp_segment_id          INTEGER NOT NULL,
                auth_seq_id             INTEGER,
                auth_seq_id_ins_code    VARCHAR,
                pdb_seq_id              INTEGER,
                unp_seq_id              INTEGER,
                observed                VARCHAR,
                dbentry_id              BIGINT,
                accession               VARCHAR,
                name                    VARCHAR,
                type                    VARCHAR,
                unp_one_letter_code     VARCHAR,
                pdb_one_letter_code     VARCHAR,
                chem_comp_id            VARCHAR,
                mh_id                   INTEGER,
                tax_id                  INTEGER,
                canonical_acc           BOOLEAN NOT NULL,
                reference_acc           VARCHAR,
                best_mapping            BOOLEAN,
                residue_id              VARCHAR
                );
            """)

    def create_sifts_xref_segment_table(self):
        """Create the ``sifts_xref_segment`` table if it does not already exist."""
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS sifts_xref_segment (
                entry_id            VARCHAR NOT NULL,
                entity_id           INTEGER NOT NULL,
                id                  INTEGER NOT NULL,
                auth_asym_id        VARCHAR,
                struct_asym_id      VARCHAR NOT NULL,
                accession           VARCHAR,
                name                VARCHAR,
                seq_version         INTEGER,
                unp_start           INTEGER,
                pdb_start           INTEGER,
                unp_end             INTEGER,
                pdb_end             INTEGER,
                auth_start          INTEGER,
                auth_start_icode    VARCHAR,
                auth_end            INTEGER,
                auth_end_icode      VARCHAR,
                conflicts           INTEGER,
                modifications       VARCHAR,
                unp_alignment       VARCHAR,
                pdb_alignment       VARCHAR,
                identity            DOUBLE,
                score               DOUBLE,
                best_mapping        BOOLEAN NOT NULL,
                canonical_acc       BOOLEAN NOT NULL,
                reference_acc       VARCHAR,
                chimera             BOOLEAN NOT NULL
            );
        """)

    def bulk_load_from_entries(self, input_dir: str) -> None:
        """Bulk-load all segment and residue CSV files under *input_dir* into DuckDB.

        The method discovers ``*_seg.csv.gz`` and ``*_res.csv.gz`` files
        recursively under *input_dir* (excluding ``*_nf90_*`` variants) and
        replaces the contents of ``sifts_xref_segment`` and
        ``sifts_xref_residue`` in one pass each.

        Args:
            input_dir: Root directory that contains the per-entry
                ``sifts/`` sub-directories produced by the segment generation
                pipeline.
        """
        base = str(Path(input_dir))
        self._bulk_load_table(
            "sifts_xref_segment", f"{base}/*_seg.csv.gz", "%_nf90_seg.csv.gz"
        )
        self._bulk_load_table(
            "sifts_xref_residue", f"{base}/*_res.csv.gz", "%_nf90_res.csv.gz"
        )
        logger.info("DuckDB bulk load complete.")

    def _bulk_load_table(
        self, table: str, glob_pattern: str, exclude: str
    ) -> None:
        """Delete and reload all rows in *table* from matching CSV files.

        Uses DuckDB's ``read_csv`` to stream compressed CSV files directly,
        inferring column names and types from the existing table schema.

        Args:
            table: Name of the target DuckDB table (e.g. ``sifts_xref_segment``).
            glob_pattern: Shell-style glob passed to DuckDB ``glob()`` to
                discover candidate files (e.g. ``/data/*_seg.csv.gz``).
            exclude: SQL ``LIKE`` pattern for files to exclude
                (e.g. ``%_nf90_seg.csv.gz``).
        """
        files = self.conn.execute(
            "SELECT list(file) FROM glob(?) WHERE NOT file LIKE ?",
            [glob_pattern, exclude],
        ).fetchone()[0]

        if not files:
            logger.warning("No files found for %s", table)
            return

        logger.info("%d files found for %s", len(files), table)

        schema = self.conn.execute(
            "SELECT column_name, data_type FROM information_schema.columns "
            "WHERE table_name = ? ORDER BY ordinal_position",
            [table],
        ).fetchall()
        col_names = [row[0] for row in schema]
        col_types = [row[1] for row in schema]

        self.conn.execute(f"DELETE FROM {table}")
        self.conn.execute(
            f"INSERT INTO {table} SELECT * FROM read_csv("
            f"?, header=false, column_names=?, column_types=?)",
            [files, col_names, col_types],
        )
        logger.info("Loaded into %s", table)
