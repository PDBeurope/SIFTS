#!/usr/bin/env python3

import argparse
from collections.abc import Mapping
from datetime import datetime
from pathlib import Path

import duckdb
import gemmi
from funcy.debug import log_durations

import pdbe_sifts.segments_generation.generate_xref_csv as generate_xref_csv
from pdbe_sifts.base.exceptions import NotAPolyPeptide, ObsoleteUniProtError
from pdbe_sifts.base.log import logger
from pdbe_sifts.base.utils import SiftsAction
from pdbe_sifts.mmcif.chem_comp import ChemCompMapping
from pdbe_sifts.mmcif.entry import Entry
from pdbe_sifts.segments_generation.alignment import helper
from pdbe_sifts.segments_generation.alignment.helper import (
    CustomSequenceAccession,
)
from pdbe_sifts.segments_generation.connectivity.process_connectivity import (
    ConnectivityCheck,
)
from pdbe_sifts.segments_generation.get_list_of_mappings import (
    get_curated_db_mappings,
)
from pdbe_sifts.unp.unp import UNP


class SiftsAlign:
    def __init__(
        self,
        cif_file,
        out_dir,
        db_conn_str=None,
        nf90_mode=False,
        unp_mode=None,
        connectivity_mode=True,
    ):
        """Initialise the per-entry SIFTS segment and residue mapping pipeline.

        Args:
            cif_file: Path to the mmCIF file to process (``.cif`` or ``.cif.gz``).
            out_dir: Root output directory; per-entry CSVs are written to
                ``{out_dir}/{entry_id}/sifts/``.
            db_conn_str: Path to a DuckDB file produced by ``global_mappings``.
                Used to retrieve the best UniProt accession per chain.
                Mutually exclusive with *unp_mode*.
            nf90_mode: When ``True``, disables the ≥ 90 % identity filter and
                keeps all alignment hits regardless of sequence identity.
            unp_mode: Manual mapping override.  Either a comma-separated
                ``auth_asym_id:accession`` string (e.g. ``"A:P00963,B:P00963"``)
                or a path to a custom FASTA file whose headers follow the
                ``>{auth_asym_id}|{sequence_id}`` convention.
                Mutually exclusive with *db_conn_str*.
            connectivity_mode: When ``True`` (default), applies the connectivity
                correction step that reassigns gap residues to adjacent segments
                when a covalent peptide bond is detected.
        """
        self.cif_file = str(cif_file)

        self.nf90_mode = nf90_mode
        self.unp_mode = unp_mode
        self.db_conn_str = db_conn_str
        self.sifts_mapping = {}
        self.out_dir = out_dir
        self.NFC = {}
        self.NFT = {}
        self.connectivity_mode = connectivity_mode
        self.custom_sequences: dict = {}
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
        # DB connection is optional: only required when no -m FASTA is given
        self.conn = (
            duckdb.connect(self.db_conn_str, read_only=True)
            if db_conn_str
            else None
        )
        logger.info("Loading chem comp three-letter to one-letter mapping")
        self.cc = ChemCompMapping()

    @log_durations(logger.debug)
    def process_entry(self, entry_id):
        cif_file = self.cif_file
        if self.no_used_cif_category_modified(cif_file):
            logger.info(
                f"{entry_id}: Modification in non-used cif category detected. Skipping."
            )
            return

        try:
            logger.info(f"Processing [{entry_id}]")
            entry = Entry(entry_id, self.cc, cif_file)
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
        chain_to_entity = {
            chain: entry.chains[chain].entity_id for chain in chain_lst
        }

        mappings = self.get_mappings(entry_id, chain_lst, chain_to_entity)
        logger.info(mappings)

        # Inject AFTER get_mappings: _parse_fasta_mapping() populates self.custom_sequences
        # inside get_mappings. Key by accession ("BLOOD"), not chain_id ("A"), because
        # get_accession(entry, acc) looks up entry.accessions[acc].
        if self.custom_sequences:
            entry.accessions.update(
                {v.accession: v for v in self.custom_sequences.values()}
            )

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
                connectivity_check = ConnectivityCheck(
                    em.chain_obj, em.repeated_acc
                )
                em.chain_obj.segments = connectivity_check.check_segments_conn()

        Path(self.out_dir).mkdir(parents=True, exist_ok=True)
        generate_xref_csv.insert_mappings(
            self.out_dir, entry, self.nf90_mode, self.conn
        )
        logger.info(f"Processed [{entry_id}]")

    def remove_existing_files(self, entry_id):
        for f in Path(self.out_dir).glob(f"{entry_id}_*.csv.gz"):
            f.unlink()

    def get_mappings(self, entry_id, chain_lst, chain_to_entity):
        mappings: Mapping[str, list[helper.SMapping]] = {}
        for chain in chain_lst:
            mappings[chain] = []
        if self.conn:
            mappings = {
                **mappings,
                **get_curated_db_mappings(
                    entry_id, chain_lst, self.conn, chain_to_entity
                ),
            }
        mappings = self._parse_user_mapping(mappings)
        return mappings

    def _parse_user_mapping(self, entry_mapping):
        if not self.unp_mode:
            return entry_mapping
        if Path(self.unp_mode).is_file():
            return self._parse_fasta_mapping(entry_mapping, Path(self.unp_mode))
        return self._parse_accession_mapping(entry_mapping)

    def _parse_accession_mapping(self, entry_mapping):
        """Legacy mode: -m 'A:P00963,B:P00963' (UniProt accessions)."""
        mapp: Mapping[str, list[helper.SMapping]] = {}
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

    def _parse_fasta_mapping(self, entry_mapping, fasta_path: Path):
        """FASTA mode: -m sequences.fasta with headers >{auth_asym_id}|{sequence_id}."""
        mapp: dict = {}

        def _flush(header: str, seq_parts: list):
            parts = header.split("|")
            chain_id = parts[0]
            accession = parts[1] if len(parts) > 1 else chain_id
            name = (
                parts[2] if len(parts) > 2 else ""
            )  # optional: >{chain}|{acc}|{name}
            sequence = "".join(seq_parts)
            self.custom_sequences[chain_id] = CustomSequenceAccession(
                accession, sequence, name
            )
            mapp.setdefault(chain_id, []).append(
                helper.SMapping(accession, 0, 0)
            )

        current_header = None
        current_seq: list[str] = []
        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header is not None:
                        _flush(current_header, current_seq)
                    current_header = line[1:]
                    current_seq = []
                elif line:
                    current_seq.append(line)
        if current_header is not None:
            _flush(current_header, current_seq)

        logger.info(f"Loaded {len(mapp)} custom sequences from {fasta_path}")
        return {**entry_mapping, **mapp}

    def _is_future_date(self, date: str) -> bool:
        return datetime.strptime(date, "%Y-%m-%d") > datetime.now()

    def no_used_cif_category_modified(self, cif_file: str) -> bool:
        if not self.used_cif_categories:
            logger.debug("No categories to check for modifications")
            return False

        self.used_cif_categories = {
            cat.lstrip("_") for cat in self.used_cif_categories
        }

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
            logger.debug(
                f"Modified categories: {', '.join(modified_categories)}"
            )
            logger.debug(
                f"Used categories: {', '.join(self.used_cif_categories)}"
            )
            logger.info("None of the modified categories are used.")
            return True

        logger.info(
            f"Modified categories found: {', '.join(modified_used_categories)}"
        )
        return False


@log_durations(logger.info)
def run():
    parser = argparse.ArgumentParser(
        "Segment generation in SIFTS, generates seg_csv, res_csv"
    )

    parser.add_argument(
        "-i",
        "--cif-input",
        required=True,
        action=SiftsAction,
        help="Input CIF file to process",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        action=SiftsAction,
        help="Output directory for CSV files",
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
            "User-defined mapping. Either UniProt accessions 'A:P00963,B:P00963', "
            "or a path to a FASTA file with headers >{auth_asym_id}|{sequence_id}."
        ),
    )

    parser.add_argument(
        "-d",
        "--duckdb",
        required=False,
        default=None,
        help=(
            "DuckDB hits file to query mappings from. Optional when -m/--mapping is provided."
        ),
    )

    parser.add_argument(
        "--entry",
        required=False,
        default=None,
        help="Entry ID to process. If omitted, derived from _entry.id in the CIF file.",
    )

    args = parser.parse_args()

    # Validate: need at least one of -d or -m
    if not args.duckdb and not args.mapping:
        parser.error("At least one of -d/--duckdb or -m/--mapping is required.")

    cif_input = args.cif_input
    if not Path(cif_input).is_file():
        parser.error(f"-i must be a CIF file, got: {cif_input}")

    if args.entry:
        entry_id = args.entry
    else:
        block = gemmi.cif.read(str(cif_input)).sole_block()
        entry_id = block.find_value("_entry.id").strip('"').lower()
        logger.info(f"Derived entry_id '{entry_id}' from mmCIF _entry.id field")

    logger.info(vars(args))
    sifts_align = SiftsAlign(
        cif_input,
        args.output_dir,
        args.duckdb,
        nf90_mode=args.nf90,
        unp_mode=args.mapping,
        connectivity_mode=args.connectivity,
    )
    sifts_align.process_entry(entry_id)
    if sifts_align.conn:
        sifts_align.conn.close()


if __name__ == "__main__":
    run()
