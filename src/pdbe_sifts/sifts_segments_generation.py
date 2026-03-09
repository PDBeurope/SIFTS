#!/usr/bin/env python3

import argparse
import duckdb
from datetime import datetime
from pathlib import Path
from collections.abc import Mapping
from funcy.debug import log_durations

import gemmi

from pdbe_sifts.base.log import logger
from pdbe_sifts.mmcif.entry import Entry
from pdbe_sifts.mmcif.chem_comp import ChemCompMapping

from pdbe_sifts.base import pdbe_path
from pdbe_sifts.mmcif.mmcif_helper import NotAPolyPeptide
from pdbe_sifts.base.batchable import Batchable
from pdbe_sifts.segments_generation.get_list_of_mappings import get_curated_db_mappings
import pdbe_sifts.segments_generation.generate_xref_csv as generate_xref_csv
from pdbe_sifts.database.sifts_db_wrapper import SiftsDB
from pdbe_sifts.segments_generation.alignment import helper
from pdbe_sifts.config import load_config

conf = load_config()


class SiftsAlign(Batchable):
    """Determine segments and related residues from a sequence alignment.

    Aligns the sequences for each chain of the entry and aligns them against
    UniProt. Identifies segments and residues and writes them to CSV files.
    Optionally writes to the DuckDB database as well.

    Note:
        When dbmode=True, DuckDB does not support concurrent writes from
        multiple processes. Use --workers 1 (or run_batch with workers=1)
        to avoid locking errors.
    """

    failure_threshold = 0.01

    def __init__(
        self,
        cif_dir,
        out_dir,
        file_duckdb,
        unp_dir,
        nf90_mode=False,
        unp_mode=None,
        dbmode=False,
    ):
        self.cif_dir = cif_dir if cif_dir else conf.location.work.data_entry_dir
        self.file_duckdb = file_duckdb
        self.unp_dir = unp_dir if unp_dir else conf.cache.uniprot
        self.nf90_mode = nf90_mode
        self.dbmode = dbmode
        self.unp_mode = unp_mode
        self.sifts_mapping = {}
        self.out_dir = out_dir
        self.NFC = {}
        self.NFT = {}
        self.entry_file_path = conf.lists.entries_all
        self.db_batch_size = 1000
        """Number of CSV files to load per DuckDB batch in teardown()."""
        self.used_cif_categories = {
            "entity_poly",
            "pdbx_struct_mod_residue",
            "pdbx_poly_seq_scheme",
            "struct_ref_seq_dif",
            "struct_ref",
            "struct_ref_seq",
            "entity_src_nat",
            "entity_src_gen",
            "pdbx_entity_src_syn",
            "entity",
            "pdbx_database_status",
            "pdbx_audit_revision_history",
        }

    def worker_setup(self):
        """Open a read-only DuckDB connection and load ChemComp mapping in each worker.

        Read-only allows multiple workers to query the DB concurrently without
        locking issues. Writes only happen in teardown() on the main process.
        """
        logger.info("Loading chem comp three-letter to one-letter mapping")
        self.cc = ChemCompMapping()
        self.conn = duckdb.connect(self.file_duckdb, read_only=True)

    def worker_teardown(self):
        self.conn.close()

    def teardown(self):
        """If dbmode is enabled, bulk-load all generated CSVs into DuckDB."""
        if not self.dbmode:
            return
        conn = duckdb.connect(self.file_duckdb)
        try:
            SiftsDB(conn).bulk_load_from_dir(self.out_dir, self.db_batch_size)
        finally:
            conn.close()

    @log_durations(logger.debug)
    def process_entry(self, entry_id):
        cif_file_path = pdbe_path.get_clean_mmcif(entry_id, base_dir=self.cif_dir)
        if self.no_used_cif_category_modified(cif_file_path):
            logger.info(
                f"{entry_id}: Modification in non-used cif category detected. Skipping."
            )
            return

        try:
            logger.info("Processing [%s]" % entry_id)
            entry = Entry(entry_id, cif_file_path)
        except NotAPolyPeptide:
            logger.warning(
                f"No pdbx_poly_seq_scheme category found for {entry_id}. Skipping"
            )
            return

        entity_lst = list(entry.entities.keys())
        if not entity_lst:
            logger.warning(f"No polypeptide chains found for entry {entry_id}")
            return
        logger.debug(f"Entities: {entity_lst}")
        mappings = self.get_mappings(entry_id, entity_lst)

        logger.info(mappings)

        for entity, entity_mapping in mappings.items():
            logger.info(f"Processing {entry_id} entity {entity}")

            em = helper.EntryMapping(
                entry, entity, entity_mapping, self.nf90_mode, self.NFT, self.NFC,
                unp_dir=self.unp_dir,
            )

            if not em.set_chain_accessions():
                logger.warning(f"Skipping {entry_id} entity {entity}")
                self.remove_existing_files(entry_id)
                continue

            em.process()

        entry_out_dir = pdbe_path.get_entry_dir(entry_id, self.out_dir, "sifts")

        Path(entry_out_dir).mkdir(parents=True, exist_ok=True)
        generate_xref_csv.insert_mappings(
            entry_out_dir, entry, self.nf90_mode, self.conn
        )
        # DuckDB writes happen in teardown() via bulk CSV load (dbmode=True).
        logger.info("Processed [%s]" % entry_id)

    def remove_existing_files(self, entry_id):
        """Remove existing CSV files for entry_id."""
        entry_out_dir = pdbe_path.get_entry_dir(entry_id, self.out_dir, "sifts")
        for f in Path(entry_out_dir).rglob("*.csv.gz"):
            f.unlink()

    def get_mappings(self, entry_id, entity_lst):
        mappings: Mapping[str, list[helper.SMapping]] = {}
        for entity in entity_lst:
            mappings[entity] = []
        mappings = {
            **mappings,
            **get_curated_db_mappings(entry_id, entity_lst, self.conn, self.unp_dir),
        }
        return mappings

    # --- CIF category change detection ---

    def no_used_cif_category_modified(self, cif_file: str) -> bool:
        """Return True if the CIF file has no modifications in the categories used.

        Scenarios:
        1. Missing revision history/category blocks → return False (process it)
        2. self.used_cif_categories is empty → return False
        3. Latest revision is not in the future → return False
        4. Modified categories exist but none are used by this task → return True (skip)
        """
        if not self.used_cif_categories:
            logger.debug("No categories to check for modifications")
            return False

        block = gemmi.cif.read(str(cif_file)).sole_block()
        history = block.find(
            "_pdbx_audit_revision_history.", ["ordinal", "revision_date"]
        )
        ordinals = [
            row["ordinal"]
            for row in history
            if self._is_future_date(row["revision_date"])
        ]

        if not ordinals:
            logger.info("No future revisions found")
            return False

        categories = block.find(
            "_pdbx_audit_revision_category.", ["revision_ordinal", "category"]
        )
        if not categories:
            logger.info("No pdbx_audit_revision_category found")
            return False

        modified_categories = {
            row["category"]
            for row in categories
            if row["revision_ordinal"] in ordinals
        }
        modified_used_categories = modified_categories.intersection(
            self.used_cif_categories
        )
        if not modified_used_categories:
            logger.debug(f"Modified categories: {', '.join(modified_categories)}")
            logger.debug(f"Used categories: {', '.join(self.used_cif_categories)}")
            logger.info("None of the modified categories are used.")
            return True

        logger.info(
            f"Modified categories found: {', '.join(modified_used_categories)}"
        )
        return False

    @staticmethod
    def _is_future_date(date: str) -> bool:
        return datetime.strptime(date, "%Y-%m-%d") > datetime.now()


@log_durations(logger.info)
def run():
    parser = argparse.ArgumentParser(
        "Segment generation in SIFTS, generates seg_csv, res_csv",
        add_help=False,
    )

    parser.add_argument(
        "-i",
        "--cif-input-dir",
        default=conf.location.work.data_entry_dir,
        help="Base location for mmCIF files",
    )
    parser.add_argument(
        "-db",
        "--db",
        required=True,
        help="duckdb file location",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=False,
        default=conf.location.work.data_entry_dir,
        help="Base location for output CSV files.",
    )
    parser.add_argument(
        "-unp",
        "--unp-dir",
        required=False,
        default=conf.cache.uniprot,
        help="Base location for unp files.",
    )
    parser.add_argument(
        "-nf90",
        "--nf90",
        action="store_true",
        default=False,
        help="UniRef90 mode (default: False)",
    )
    parser.add_argument(
        "-w",
        "--write-to-db",
        action="store_true",
        default=False,
        help="Additionally write to duckdb file (default: False)",
    )
    parser.add_argument(
        "-m",
        "--mapping",
        help=(
            "User defined uniprot accession given (default: False), example use: "
            "python sifts_segments_generation.py -m A:P00963,B:P00963 ..."
        ),
    )

    # Parse only the custom arguments; let Batchable.main() handle single/batch
    custom_args, remaining = parser.parse_known_args()

    sifts_align = SiftsAlign(
        custom_args.cif_input_dir,
        custom_args.output_dir,
        custom_args.db,
        custom_args.unp_dir,
        nf90_mode=custom_args.nf90,
        unp_mode=custom_args.mapping,
        dbmode=custom_args.write_to_db,
    )
    sifts_align.main(remaining)


if __name__ == "__main__":
    run()
