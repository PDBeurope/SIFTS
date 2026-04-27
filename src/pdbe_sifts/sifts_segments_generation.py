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
from pdbe_sifts.base.utils import SiftsAction, validate_entry_id
from pdbe_sifts.mmcif.chem_comp import ChemCompMapping
from pdbe_sifts.mmcif.curator.cif_sequence_update import CifSequenceUpdater
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
    """Orchestrate the SIFTS segment and residue mapping pipeline for one PDB entry.

    Reads an mmCIF file, resolves UniProt accessions for every polypeptide
    chain (from a DuckDB hits database and/or a user-supplied mapping),
    runs pairwise sequence alignments, applies optional connectivity
    correction, and writes gzipped segment and residue CSV files to the
    specified output directory.

    The class is designed to be instantiated once per batch and reused across
    multiple calls to :meth:`process_entry`.

    Attributes:
        cif_file: Absolute path to the mmCIF file being processed.
        nf90_mode: ``True`` when the ≥ 90 % identity filter is disabled.
        unp_mode: Raw value of the ``-m`` / ``--mapping`` argument, or
            ``None`` if not supplied.
        db_conn_str: Path to the DuckDB file, or ``None``.
        sifts_mapping: Internal mapping state (populated during processing).
        out_dir: Root output directory for CSV files.
        NFC: Non-fragment connectivity cache.
        NFT: Non-fragment topology cache.
        connectivity_mode: ``True`` when connectivity correction is enabled.
        custom_sequences: Dict of chain-ID → :class:`CustomSequenceAccession`
            populated from a user-supplied FASTA file.
        used_cif_categories: Set of mmCIF category names monitored for
            future-dated revisions.
        conn: Active ``duckdb.DuckDBPyConnection``, or ``None``.
        cc: :class:`~pdbe_sifts.mmcif.chem_comp.ChemCompMapping` instance.
    """

    def __init__(
        self,
        cif_file,
        out_dir,
        db_conn_str=None,
        nf90_mode=False,
        unp_mode=None,
        connectivity_mode=True,
        tax_tsv=None,
    ):
        """Initialise the per-entry SIFTS segment and residue mapping pipeline.

        Args:
            cif_file: Path to the mmCIF file to process (``.cif`` or ``.cif.gz``).
            out_dir: Root output directory; per-entry CSVs are written to
                ``{out_dir}/{entry_id}/sifts/``.
            db_conn_str: Path to a DuckDB file produced by ``sequence_match``.
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
            tax_tsv: Optional path to a TSV file overriding the taxonomy ID per
                entity (columns: ``entity_id``, ``tax_id``).  When supplied,
                each chain whose entity_id appears in the TSV will have its
                ``tax_id`` replaced before alignment.
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

        # Build entity_id → tax_id override map from TSV
        self.tax_override: dict[str, int] = {}
        if tax_tsv:
            import pandas as _pd

            df = _pd.read_csv(
                tax_tsv, sep="\t", dtype={"entity_id": str, "tax_id": int}
            )
            self.tax_override = dict(
                zip(df["entity_id"], df["tax_id"], strict=False)
            )
            logger.info(
                "Loaded %d tax_id overrides from %s",
                len(self.tax_override),
                tax_tsv,
            )

        logger.info("Loading chem comp three-letter to one-letter mapping")
        self.cc = ChemCompMapping()

    @log_durations(logger.debug)
    def _ensure_poly_seq_scheme(self, cif_path: str) -> str:
        """Return *cif_path* unchanged if ``_pdbx_poly_seq_scheme`` is present.

        If the category is missing, generates an enriched CIF via
        :class:`~pdbe_sifts.mmcif.curator.cif_sequence_update.CifSequenceUpdater`
        and returns the new path.  The enriched file is written alongside the
        input with a ``_pdbx_added`` suffix, e.g.::

            9b4h_updated.cif    → 9b4h_updated_pdbx_added.cif
            9b4h_updated.cif.gz → 9b4h_updated_pdbx_added.cif
        """
        block = gemmi.cif.read(str(cif_path)).sole_block()
        if block.get_mmcif_category("_pdbx_poly_seq_scheme"):
            return cif_path

        p = Path(cif_path)
        if p.name.endswith(".cif.gz"):
            enriched = Path(p.parent, (p.name[:-7] + "_pdbx_added.cif"))
        else:
            enriched = p.with_name(p.stem + "_pdbx_added" + p.suffix)

        logger.info(
            "_pdbx_poly_seq_scheme missing in %s — generating enriched CIF: %s",
            p.name,
            enriched.name,
        )
        CifSequenceUpdater(str(cif_path), str(enriched)).process()
        return str(enriched)

    def process_entry(self, entry_id):
        """Run the full SIFTS segment and residue mapping pipeline for one entry.

        Parses the mmCIF file, resolves UniProt mappings for each chain,
        runs pairwise alignments, applies connectivity correction, and
        writes output CSV files to ``self.out_dir``.

        Args:

            entry_id: PDB identifier (e.g. ``"1abc"``).

        Raises:
            ValueError: If *entry_id* fails the path-safety allowlist check.
        """
        validate_entry_id(entry_id)
        cif_file = self.cif_file
        if self.no_used_cif_category_modified(cif_file):
            logger.info(
                f"{entry_id}: Modification in non-used cif category detected. Skipping."
            )
            return

        cif_file = self._ensure_poly_seq_scheme(cif_file)

        try:
            logger.info(f"Processing [{entry_id}]")
            entry = Entry(entry_id, self.cc, cif_file)
        except NotAPolyPeptide:
            logger.warning(
                f"No pdbx_poly_seq_scheme category found for {entry_id}. Skipping"
            )
            return

        if self.tax_override:
            for chain in entry.chains.values():
                if chain.entity_id in self.tax_override:
                    chain.tax_id = self.tax_override[chain.entity_id]
                    logger.debug(
                        "Overriding tax_id for entity %s chain %s → %d",
                        chain.entity_id,
                        chain.auth_asym_id,
                        chain.tax_id,
                    )

        chain_lst = list(entry.chains.keys())
        if not chain_lst:
            logger.warning(f"No polypeptide chains found for entry {entry_id}")
            return
        logger.debug(f"Chains: {chain_lst}")
        chain_to_entity = {
            chain: entry.chains[chain].entity_id for chain in chain_lst
        }

        mappings = self.get_mappings(entry_id, chain_lst, chain_to_entity)
        logger.debug(mappings)

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
        """Delete any existing output CSVs for *entry_id* in the output directory.

        Args:
            entry_id: PDB entry identifier; used to glob ``{entry_id}_*.csv.gz``.
        """
        for f in Path(self.out_dir).glob(f"{entry_id}_*.csv.gz"):
            f.unlink()

    def get_mappings(self, entry_id, chain_lst, chain_to_entity):
        """Retrieve UniProt mappings for all chains of an entry.

        Initialises every chain with an empty mapping list, then populates
        it from the DuckDB hits database (if available) and overlays any
        user-supplied mapping from ``self.unp_mode``.

        Args:
            entry_id: PDB entry identifier.
            chain_lst: List of author chain IDs present in the entry.
            chain_to_entity: Mapping of chain ID → mmCIF entity ID.

        Returns:
            Dict mapping chain ID → list of :class:`~helper.SMapping`.
        """
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
        """Overlay user-supplied mappings onto the base mapping dict.

        Dispatches to ``_parse_fasta_mapping`` when ``self.unp_mode`` points
        to a file, or to ``_parse_accession_mapping`` for a comma-separated
        ``chain:accession`` string.  Returns *entry_mapping* unchanged when
        ``self.unp_mode`` is ``None``.

        Args:
            entry_mapping: Base chain→SMapping dict (modified in-place via
                dict merge).

        Returns:
            Updated mapping dict with user overrides applied.
        """
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
        """FASTA mode: -m sequences.fasta with headers >{entry_id}|{auth_asym_id}|{sequence_id}."""
        mapp: dict = {}

        def _flush(header: str, seq_parts: list):
            """Finalise one FASTA record and register it in *mapp*.

            Called each time a new ``>`` header line is encountered (to emit
            the preceding record) and once more after the last line of the
            file.  Populates ``self.custom_sequences`` and appends an
            :class:`~helper.SMapping` entry to *mapp*.

            Args:
                header: FASTA header text (without the leading ``>``), in the
                    format ``{entry_id}|{auth_asym_id}|{sequence_id}`` or
                    ``{entry_id}|{auth_asym_id}|{sequence_id}|{name}``.
                seq_parts: List of non-empty sequence lines accumulated since
                    the last header; joined to form the full sequence string.
            """
            parts = header.split("|")
            # parts[0] = entry_id (ignored in single-entry mode — entry already known)
            chain_id = parts[1] if len(parts) > 1 else ""
            accession = parts[2] if len(parts) > 2 else chain_id
            name = parts[3] if len(parts) > 3 else ""
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
        """Return ``True`` if *date* is strictly after today.

        Args:
            date: Date string in ``YYYY-MM-DD`` format.
        """
        return datetime.strptime(date, "%Y-%m-%d") > datetime.now()

    def no_used_cif_category_modified(self, cif_file: str) -> bool:
        """Check whether any SIFTS-relevant mmCIF category has a future revision.

        Reads ``_pdbx_audit_revision_history`` and
        ``_pdbx_audit_revision_category`` to find revisions dated in the
        future that touch categories used by the pipeline.

        Args:
            cif_file: Path to the mmCIF file to inspect.

        Returns:
            ``True`` if a future revision modifies a used category (the
            entry should be skipped until the revision date is reached).
            ``False`` if no such revision exists or the category list is
            empty.
        """
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
    """CLI entry point for the SIFTS segment generation pipeline.

    Parses command-line arguments, constructs a :class:`SiftsAlign` instance,
    and calls :meth:`SiftsAlign.process_entry` for the target entry.  The
    entry ID is either taken from the ``--entry`` argument or derived from
    the ``_entry.id`` field of the mmCIF file.

    At least one of ``-d`` / ``--duckdb`` or ``-m`` / ``--mapping`` must be
    provided; the function calls ``parser.error`` and exits otherwise.

    Command-line arguments:
        -i / --cif-input: Path to the input mmCIF file.
        -o / --output-dir: Root output directory for CSV files.
        -nf90 / --nf90: Enable NF90 mode (default: ``False``).
        --no-connectivity: Disable connectivity correction (default: enabled).
        -m / --mapping: User-defined chain→accession mapping or FASTA file.
        -d / --duckdb: Path to the DuckDB hits file.
        --entry: PDB entry ID override; derived from the mmCIF file if omitted.
    """
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
            "or a path to a FASTA file with headers >{entry_id}|{auth_asym_id}|{sequence_id}."
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
        entry_id = validate_entry_id(args.entry.lower())
    else:
        block = gemmi.cif.read(str(cif_input)).sole_block()
        raw_id = block.find_value("_entry.id").strip('"').lower()
        entry_id = validate_entry_id(raw_id)
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
