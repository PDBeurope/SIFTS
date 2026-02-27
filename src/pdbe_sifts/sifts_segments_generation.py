#!/usr/bin/env python3

import argparse
import io
import csv
import json
import duckdb
from timeit import default_timer as timer
from pathlib import Path
import pickle
from collections.abc import Mapping
from funcy.debug import log_durations

from pdbe_sifts.base.log import logger
from pdbe_sifts.mmcif.entry import Entry
from pdbe_sifts.base.utils import get_date, make_path
from pdbe_sifts.mmcif.chem_comp import ChemCompMapping

from pdbe_sifts.base.parser import parse_with_base_parser
from pdbe_sifts.base import pdbe_path
from pdbe_sifts.mmcif.mmcif_helper import NotAPolyPeptide
from pdbe_sifts.base.batchable import Batchable
from pdbe_sifts.segments_generation.get_list_of_mappings import get_curated_db_mappings
import pdbe_sifts.segments_generation.generate_xref_csv as generate_xref_csv
from pdbe_sifts.database.sifts_db_wrapper import SiftsDB
from pdbe_sifts.base.exceptions import ObsoleteUniProtError
from pdbe_sifts.segments_generation.alignment import helper
from pdbe_sifts.unp.unp import UNP
from pdbe_sifts.config import load_config
# from orc.sifts.uniref90_pkl import NF90Coverage, NF90TaxID

conf = load_config()

class SiftsAlign(Batchable):
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
        """Determine segments and related residues from a sequence alignment.

        Aligns the sequences for each chain of the entry and aligns them against uniprot.
        Also identifies segments and residues and writes them to CSV files
        Optionally writes to the database as well.

        Args:
            nf90_mode(bool): Whether processing NF90 or normal ISOFORM SIFTS
            unp_mode(str): Provide mapping A:P00963,B:P00963
        """

        self.cif_dir = cif_dir if cif_dir else conf.location.work.data_entry_dir
        self.file_duckdb = file_duckdb
        self.conn = duckdb.connect(self.file_duckdb)
        self.unp_dir = unp_dir if unp_dir else conf.cache.uniprot
        self.nf90_mode = nf90_mode
        self.dbmode = dbmode
        self.failure_threshold = 0.01
        self.unp_mode = unp_mode
        self.sifts_mapping = {}
        self.out_dir = out_dir
        self.NFC = {}
        self.NFT = {}
        self.initial_memory = 8000
        self.retry_memory = 16000
        self.entry_file_path = conf.lists.entries_all
        self.used_cif_categories = [
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
        ]

    def before_job_start(self):
        logger.info("Loading chem comp three-letter to one-letter mapping")
        self.cc = ChemCompMapping()
    
    def after_job_end(self):
        self.conn.close()

    @log_durations(logger.debug)
    def process_entry(self, entry_id):
        cif_file_path = pdbe_path.get_clean_mmcif(entry_id, base_dir=self.cif_dir)
        if self.no_used_cif_category_modified(
            cif_file_path
        ):
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
                entry, entity, entity_mapping, self.nf90_mode, self.NFT, self.NFC
            )

            if not em.set_chain_accessions():
                logger.warning(f"Skipping {entry_id} chain {chain}")
                self.remove_existing_files(entry_id)
                continue

            em.process()

        #     if self.nf90_mode:
        #         try:
        #             entry.chains[chain].skip = True
        #         except KeyError:
        #             # the chain was not found in the mmcif or it is not
        #             # a polypeptide. It is not going into the DB anyway
        #             pass

        entry_out_dir = pdbe_path.get_entry_dir(entry_id, self.out_dir, "sifts")

        Path(entry_out_dir).mkdir(parents=True, exist_ok=True)
        segments, residues = generate_xref_csv.insert_mappings(
            entry_out_dir, entry, self.nf90_mode, self.conn
        )
        if self.dbmode:
            logger.info(f'Inserting mappings into {self.file_duckdb}.')
            db_wrapper = SiftsDB(self.conn)
            db_wrapper.insert_xref_segments(segments)
            db_wrapper.insert_xref_residues(residues)
            logger.info(f'Mappings inserted into {self.file_duckdb}.')

        logger.info("Processed [%s]" % entry_id)


    def remove_existing_files(self, entry_id):
        """Remove existing files for entry_id"""
        entry_out_dir = pdbe_path.get_entry_dir(entry_id, self.out_dir, "sifts")
        for f in Path(entry_out_dir).rglob("*.csv.gz"):
            f.unlink()

    def get_mappings(self, entry_id, entity_lst):
        mappings: Mapping[str, list[helper.SMapping]] = {}
        # Initialize all entities
        for entity in entity_lst:
            mappings[entity] = []
        # Replace with database mappings
        mappings = {
            **mappings,
            **get_curated_db_mappings(entry_id, entity_lst, self.conn, self.unp_dir),
        }
        # Overwrite with user specified mappins
        # mappings = self._parse_user_mapping(mappings)
        return mappings

    # def _parse_user_mapping(self, entry_mapping):
    #     mapp: Mapping[str, list[helper.SMapping]] = {}
    #     if self.unp_mode:
    #         chains = self.unp_mode.split(",")
    #         for chain in chains:
    #             chain, acc = chain.split(":")
    #             try:
    #                 unp = UNP(acc)
    #                 mapp.setdefault(chain, []).append(
    #                     helper.SMapping(unp.accession, 0, 0)
    #                 )
    #             except ObsoleteUniProtError:
    #                 logger.warning(
    #                     f"Obsolete UniProt accession provided by user: {acc}. Will be ignored"
    #                 )
    #                 continue

    #     return {**entry_mapping, **mapp}



@log_durations(logger.info)
def run():
    parser = argparse.ArgumentParser(
        "Segment generation in SIFTS, generates seg_csv, res_csv"
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
            "User defined uniprot accession given (default: False), example use:"
            "python test_sifts_alignments.py -d -wu  A:P00963,B:P00963"
        ),
    )

    # adding single,batch options
    args = parse_with_base_parser(parser)

    logger.info(vars(args))
    sifts_align = SiftsAlign(
        args.cif_input_dir,
        args.output_dir,
        args.db,
        args.unp_dir,
        nf90_mode=args.nf90,
        unp_mode=args.mapping,
        dbmode=args.write_to_db,
    )
    sifts_align.before_job_start()
    sifts_align.process_entry(args.entry)
    sifts_align.after_job_end()


if __name__ == "__main__":
    run()