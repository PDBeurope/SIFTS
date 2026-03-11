#!/usr/bin/env python3

import argparse
import duckdb
from collections.abc import Mapping
from datetime import datetime
from pathlib import Path

import gemmi
from funcy.debug import log_durations

from pdbe_sifts.base import pdbe_path
from pdbe_sifts.base.exceptions import ObsoleteUniProtError
from pdbe_sifts.base.log import logger
from pdbe_sifts.base.utils import SiftsAction
from pdbe_sifts.config import load_config
from pdbe_sifts.mmcif.chem_comp import ChemCompMapping
from pdbe_sifts.mmcif.entry import Entry
from pdbe_sifts.mmcif.mmcif_helper import NotAPolyPeptide
import pdbe_sifts.segments_generation.generate_xref_csv as generate_xref_csv
from pdbe_sifts.segments_generation.alignment import helper
from pdbe_sifts.segments_generation.connectivity.process_connectivity import ConnectivityCheck
from pdbe_sifts.segments_generation.get_list_of_mappings import get_curated_db_mappings
from pdbe_sifts.unp.unp import UNP

conf = load_config()


class SiftsAlign:
    def __init__(
        self,
        cif_dir,
        out_dir,
        db_conn_str,
        nf90_mode=False,
        unp_mode=None,
        connectivity_mode=True,
    ):
        self.cif_dir = cif_dir
        self.input = cif_dir
        self.nf90_mode = nf90_mode
        self.unp_mode = unp_mode
        self.db_conn_str = db_conn_str
        self.sifts_mapping = {}
        self.out_dir = out_dir
        self.NFC = {}
        self.NFT = {}
        self.connectivity_mode = connectivity_mode
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
        self.conn = duckdb.connect(self.db_conn_str, read_only=True)
        logger.info("Loading chem comp three-letter to one-letter mapping")
        self.cc = ChemCompMapping()

    @log_durations(logger.debug)
    def process_entry(self, entry_id):
        if self.no_used_cif_category_modified(
            pdbe_path.get_clean_mmcif(entry_id, base_dir=self.input)
        ):
            logger.info(
                f"{entry_id}: Modification in non-used cif category detected. Skipping."
            )
            return

        try:
            logger.info("Processing [%s]" % entry_id)
            entry = Entry(entry_id, self.cc, self.cif_dir)
        except NotAPolyPeptide:
            logger.warning(
                f"No pdbx_poly_seq_scheme category found for {entry_id}. Skipping"
            )
            return

        chain_lst = list(entry.chains.keys())
        if not chain_lst:
            logger.warning(f"No polypeptide chains found for entry {entry_id}")
            return
        logger.debug(f"Chains: {chain_lst}")
        chain_to_entity = {chain: entry.chains[chain].entity_id for chain in chain_lst}
        mappings = self.get_mappings(entry_id, chain_lst, chain_to_entity)
        logger.info(mappings)

        for chain, chain_mapping in mappings.items():
            logger.info(f"Processing {entry_id} chain {chain}")

            em = helper.EntryMapping(
                entry,
                chain,
                chain_mapping,
                self.nf90_mode,
                self.NFT,
                self.NFC,
                self.connectivity_mode,
            )

            if not em.set_chain_accessions():
                logger.warning(f"Skipping {entry_id} chain {chain}")
                self.remove_existing_files(entry_id)
                continue

            em.process()
            if not em.chain_obj.is_chimera and self.connectivity_mode:
                connectivity_check = ConnectivityCheck(em.chain_obj, em.repeated_acc)
                em.chain_obj.segments = connectivity_check.check_segments_conn()

        entry_out_dir = pdbe_path.get_entry_dir(entry_id, self.out_dir, "sifts")

        Path(entry_out_dir).mkdir(parents=True, exist_ok=True)
        generate_xref_csv.insert_mappings(
            entry_out_dir, entry, self.nf90_mode, self.conn
        )
        logger.info("Processed [%s]" % entry_id)

    def remove_existing_files(self, entry_id):
        entry_out_dir = pdbe_path.get_entry_dir(entry_id, self.out_dir, "sifts")
        for f in Path(entry_out_dir).rglob("*.csv.gz"):
            f.unlink()

    def get_mappings(self, entry_id, chain_lst, chain_to_entity):
        mappings: Mapping[str, list[helper.SMapping]] = {}
        for chain in chain_lst:
            mappings[chain] = []
        mappings = {
            **mappings,
            **get_curated_db_mappings(entry_id, chain_lst, self.conn, chain_to_entity),
        }
        mappings = self._parse_user_mapping(mappings)
        return mappings

    def _parse_user_mapping(self, entry_mapping):
        mapp: Mapping[str, list[helper.SMapping]] = {}
        if self.unp_mode:
            chains = self.unp_mode.split(",")
            for chain in chains:
                chain, acc = chain.split(":")
                try:
                    unp = UNP(acc)
                    mapp.setdefault(chain, []).append(
                        helper.SMapping(unp.accession, 0, 0)
                    )
                except ObsoleteUniProtError:
                    logger.warning(
                        f"Obsolete UniProt accession provided by user: {acc}. Will be ignored"
                    )
                    continue

        return {**entry_mapping, **mapp}

    def _is_future_date(self, date: str) -> bool:
        return datetime.strptime(date, "%Y-%m-%d") > datetime.now()

    def no_used_cif_category_modified(self, cif_file: str) -> bool:
        if not self.used_cif_categories:
            logger.debug("No categories to check for modifications")
            return False

        self.used_cif_categories = {cat.lstrip("_") for cat in self.used_cif_categories}

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
            row["category"] for row in categories if row["revision_ordinal"] in ordinals
        }
        modified_used_categories = modified_categories.intersection(
            self.used_cif_categories
        )
        if not modified_used_categories:
            logger.debug(f"Modified categories: {', '.join(modified_categories)}")
            logger.debug(f"Used categories: {', '.join(self.used_cif_categories)}")
            logger.info("None of the modified categories are used.")
            return True

        logger.info(f"Modified categories found: {', '.join(modified_used_categories)}")
        return False


@log_durations(logger.info)
def run():
    parser = argparse.ArgumentParser(
        "Segment generation in SIFTS, generates seg_csv, res_csv"
    )

    parser.add_argument(
        "-i",
        "--cif-input-dir",
        required=True,
        action=SiftsAction,
        default=conf.location.work.data_entry_dir,
        help="Base location for mmCIF files",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        action=SiftsAction,
        default=conf.location.work.data_entry_dir,
        help="Base location for output CSV files.",
    )

    parser.add_argument(
        "-nf90",
        "--nf90",
        action="store_true",
        default=False,
        help="UniRef90 mode (default: False)",
    )

    parser.add_argument(
        "--no-connectivity",
        dest="connectivity",
        action="store_false",
        default=True,
        help="disable connectivity mode (default: enabled)",
    )

    parser.add_argument(
        "-m",
        "--mapping",
        help=(
            "User defined uniprot accession given (default: False), example use:"
            "python test_sifts_alignments.py -d -wu  A:P00963,B:P00963"
        ),
    )

    parser.add_argument(
        "-d",
        "--duckdb",
        required=True,
        help=(
            "Database used to query mappings in single mode "
            "(and write in db tables if -w supplied)"
        ),
    )

    parser.add_argument(
        "--entry",
        required=True,
        help="Entry ID to process.",
    )

    args = parser.parse_args()

    logger.info(vars(args))
    sifts_align = SiftsAlign(
        args.cif_input_dir,
        args.output_dir,
        args.duckdb,
        nf90_mode=args.nf90,
        unp_mode=args.mapping,
        connectivity_mode=args.connectivity,
    )
    sifts_align.process_entry(args.entry)
    sifts_align.conn.close()


if __name__ == "__main__":
    run()
