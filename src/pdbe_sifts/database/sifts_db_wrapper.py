from pathlib import Path
import duckdb
import pandas as pd

from pdbe_sifts.base.log import logger

def segment_to_dict(seg):
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
    def __init__(self, conn):
        self.conn = conn
        self.create_sifts_xref_residue_table()
        self.create_sifts_xref_segment_table()


    def create_sifts_xref_residue_table(self):
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
                dbentry_id              INTEGER,
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


    def clean_segments(self, entry_ids):
        entry_id = [entry.lower() for entry in entry_ids]
        self.conn.execute(
            """
            DELETE FROM sifts_xref_segment
            WHERE entry_id IN ?
            """,
            [entry_ids]
        )


    def clean_residues(self, entry_ids):
        entry_id = [entry.lower() for entry in entry_ids]
        self.conn.execute(
            """
            DELETE FROM sifts_xref_residue 
            WHERE entry_id IN ?
            """,
            [entry_ids]
        )

    def insert_xref_segments(self, segments):
        if not segments:
            return

        rows = [segment_to_dict(s) for s in segments]
        df = pd.DataFrame(rows)
        entries = df['entry_id'].to_list()
        self.clean_segments(entries)


        self.conn.register("tmp_segments", df)
        self.conn.execute("""
            INSERT INTO sifts_xref_segment
            SELECT * FROM tmp_segments
        """)
        self.conn.unregister("tmp_segments")
    
    def insert_xref_residues(self, residues):
        if not residues:
            return

        rows = [residue_to_dict(r) for r_list in residues for r in r_list]
        df = pd.DataFrame(rows)
        entries = df['entry_id'].to_list()
        self.clean_residues(entries)

        self.conn.register("tmp_residues", df)
        self.conn.execute("""
            INSERT INTO sifts_xref_residue
            SELECT * FROM tmp_residues
        """)
        self.conn.unregister("tmp_residues")
    
    def clean_by_entry(self, entry):
        pass

    def bulk_load_from_dir(self, out_dir, batch_size: int = 1000) -> None:
        """Bulk-load all segment and residue CSV files from out_dir into DuckDB.

        Scans out_dir recursively for *_seg.csv.gz and *_res.csv.gz files,
        then loads them in batches. Each batch is idempotent: existing rows
        for those entry_ids are deleted before inserting the new ones.

        Args:
            out_dir: Root directory to scan for CSV files.
            batch_size: Number of CSV files to load per DuckDB transaction.
        """
        out_path = Path(out_dir)
        seg_files = sorted(out_path.rglob("*_seg.csv.gz"))
        res_files = sorted(out_path.rglob("*_res.csv.gz"))

        if not seg_files and not res_files:
            logger.warning("bulk_load_from_dir: no CSV files found in %s", out_dir)
            return

        logger.info(
            "Bulk loading %d segment files and %d residue files from %s",
            len(seg_files), len(res_files), out_dir,
        )
        self._bulk_load_table("sifts_xref_segment", seg_files, batch_size)
        self._bulk_load_table("sifts_xref_residue", res_files, batch_size)
        logger.info("DuckDB bulk load complete.")

    def _bulk_load_table(self, table: str, files: list, batch_size: int) -> None:
        """Load a list of gzipped CSV files into a DuckDB table in batches.

        For each batch:
          - entry_ids are extracted from filenames (e.g. 1abc_seg.csv.gz → 1abc)
          - existing rows for those entry_ids are deleted first (idempotent)
          - new rows are inserted via DuckDB's read_csv_auto()

        Args:
            table: Target table name in DuckDB.
            files: List of Path objects pointing to gzipped CSV files.
            batch_size: Number of files per batch.
        """
        total = len(files)
        if not total:
            return

        logger.info("Loading %d files into %s...", total, table)

        # Fetch the ordered column names from the table schema so that
        # read_csv_auto does not try to infer them from the first data row.
        schema = self.conn.execute(
            "SELECT column_name FROM information_schema.columns "
            "WHERE table_name = ? ORDER BY ordinal_position",
            [table],
        ).fetchall()
        col_names_sql = "[" + ", ".join(f"'{row[0]}'" for row in schema) + "]"

        for i in range(0, total, batch_size):
            batch = files[i : i + batch_size]

            entry_ids = [f.name.split("_")[0] for f in batch]
            entry_ids_sql = ", ".join(f"'{e}'" for e in entry_ids)
            file_list_sql = "[" + ", ".join(f"'{f}'" for f in batch) + "]"

            self.conn.execute(
                f"DELETE FROM {table} WHERE entry_id IN ({entry_ids_sql})"
            )
            self.conn.execute(
                f"INSERT INTO {table} SELECT * FROM read_csv_auto("
                f"{file_list_sql}, header=false, column_names={col_names_sql})"
            )

            loaded = min(i + batch_size, total)
            logger.info("  %s: %d/%d files loaded", table, loaded, total)
