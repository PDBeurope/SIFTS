#!/usr/bin/env python3
"""
SIFTS global mappings generator.

This module performs global sequence–structure mappings for PDB entries
using MMseqs2 or BLASTP. It extracts sequences from mmCIF files, builds
temporary FASTA files, runs alignments, and parses the resulting hits
into ranked mappings.
"""

import argparse
from timeit import default_timer as timer
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

from pdbe_sifts.base.log import logger
from pdbe_sifts.mmcif.entry import Entry
from pdbe_sifts.global_mappings.mmseqs_search import MmSearch
from pdbe_sifts.global_mappings.blastp import BlastP
from pdbe_sifts.global_mappings.global_mappings_parser import GlobMappingsParser
from pdbe_sifts.base.utils import get_date, make_path


class SiftsGlobalMappings:
    """Main class for generating structure–sequence mappings."""

    def __init__(
        self,
        cif_file: str | Path,
        out_dir: str | Path,
        db_file: str | Path,
        unp_csv: Optional[str | Path],
        tool: str = "mmseqs",
        threads: int = 1,
        batch_size: int = 100000,
    ):
        """
        Initialize SiftsGlobalMappings.

        Args:
            cif_file: Path to a mmCIF file or a text file listing mmCIF paths.
            out_dir: Directory where results will be written.
            db_file: Path to the preformatted sequence database (MMseqs or BLAST).
            tool: Alignment tool ('mmseqs' or 'blastp'). Default: 'mmseqs'.
            threads: Number of CPU threads to use.
            batch_size: Number of CIF files per batch when processing a .txt list.
        """
        self.cif_file = Path(cif_file)
        self.out_dir = Path(out_dir)
        self.db_file = Path(db_file)
        self.unp_csv = Path(unp_csv) if unp_csv is not None else None
        self.tool = tool
        self.threads = threads
        self.batch_size = batch_size
        self.date = get_date()

        # Output directories
        self.fasta_dir = self.out_dir / "fasta_files"
        self.fasta_dir.mkdir(parents=True, exist_ok=True)

        self.result_file_path = {}


    def generate_fasta(
        self,
        entity_seq_tax_dict: dict,
        entry_name: str,
        file_name: str,
        mode: str = "w",
    ) -> Path:
        """
        Write a temporary FASTA file for an entry or batch of entries.

        Args:
            entity_seq_tax_dict: Mapping of entity_id → (sequence, taxid).
            entry_name: Entry name (e.g. '8kd1').
            file_name: Base name for the output FASTA file.
            mode: File open mode ('w' or 'a').

        Returns:
            Path to the generated FASTA file.
        """
        fasta_path = self.fasta_dir / f"tmp_{file_name}.fasta"
        with open(fasta_path, mode) as fh:
            for entity, (seq, tax_id) in entity_seq_tax_dict.items():
                header = f">pdb|{entry_name}-{entity}|OX={tax_id}"
                fh.write(f"{header}\n{seq}\n")
        return fasta_path


    def _process_single_file(self, file_path: str | Path):
        """Extract sequence/taxonomy info from one mmCIF file."""
        file_path = Path(file_path)
        if not file_path.exists():
            return None, f"Missing CIF file: {file_path}"

        entry_name = file_path.stem
        try:
            entry = Entry(entry_name, str(file_path))
            entity_seq_tax = entry.get_entity_seq_tax()
            return (entry_name, entity_seq_tax), None
        except Exception as e:
            return None, f"Failed to process {file_path}: {e}"

    def process_input_file(self) -> tuple[Path, str]:
        """
        Process a mmCIF or text file input.

        Returns:
            (Path to temporary FASTA file, entry name or list name)
        """
        if not self.cif_file.exists():
            raise FileNotFoundError(f"Input file not found: {self.cif_file}")

        ext = self.cif_file.suffix.lower()

        # Single CIF file
        if ext == ".cif":
            entry_name = self.cif_file.stem
            entry = Entry(entry_name, str(self.cif_file))
            entity_seq_tax = entry.get_entity_seq_tax()
            fasta_path = self.generate_fasta(entity_seq_tax, entry_name, entry_name)
            return fasta_path, entry_name

        # List of CIF files
        if ext == ".txt":
            list_name = self.cif_file.stem
            fasta_path = self.fasta_dir / f"tmp_{list_name}.fasta"
            fasta_path.unlink(missing_ok=True)

            file_paths = [
                line.strip()
                for line in self.cif_file.read_text().splitlines()
                if line.strip()
            ]

            for i in range(0, len(file_paths), self.batch_size):
                batch = file_paths[i : i + self.batch_size]
                results = []

                with ProcessPoolExecutor(max_workers=self.threads) as pool:
                    futures = {pool.submit(self._process_single_file, fp): fp for fp in batch}
                    for fut in as_completed(futures):
                        res, warn = fut.result()
                        if warn:
                            logger.warning(warn)
                            continue
                        if res:
                            entry_name, entity_seq_tax = res
                            results.append((entry_name, entity_seq_tax))

                for entry_name, entity_seq_tax in results:
                    fasta_path = self.generate_fasta(entity_seq_tax, entry_name, list_name, mode="a")

            return fasta_path, list_name

        raise ValueError(f"Unsupported input format: {ext} (must be .cif or .txt)")


    def mmseqs_search(self, entry_id: str, fasta_path: Path):
        """Run an MMseqs2 search."""
        output_path = make_path(self.out_dir, entry_id, "mmseqs", f"hits_{entry_id}.tsv", self.date)
        tmp_dir = output_path.parent / f"tmp_{entry_id}"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        search = MmSearch(fasta_path, self.db_file, output_path, tmp_dir, self.threads, db_load_mode=2)
        search.run()
        self.result_file_path[entry_id] = output_path

    def blastp_search(self, entry_id: str, fasta_path: Path):
        """Run a BLASTP search."""
        output_path = make_path(self.out_dir, entry_id, "blastp", f"hits_{entry_id}.tsv", self.date)
        search = BlastP(fasta_path, self.db_file, output_path, threads=self.threads)
        search.run()
        self.result_file_path[entry_id] = output_path

    def search(self, entry_id: str, fasta_path: Path):
        """Dispatch to the selected search tool."""
        if self.tool == "mmseqs":
            self.mmseqs_search(entry_id, fasta_path)
        elif self.tool == "blastp":
            self.blastp_search(entry_id, fasta_path)
        else:
            raise ValueError(f"Unsupported tool: {self.tool}")


    def process(self):
        """Run the complete mapping pipeline."""
        logger.info(f"Processing [{self.cif_file}]")
        start = timer()
        fasta_path, entry_name = self.process_input_file()
        self.search(entry_name, fasta_path)
        logger.info(f"Parsing {self.tool} hits.")
        self.mappings = GlobMappingsParser(
            self.tool,
            self.result_file_path[entry_name],
            Path(self.result_file_path[entry_name]).parent,
            unp_csv=self.unp_csv,
        ).parse()

        end = timer()
        logger.info(f"Total time (from parsing to ranked mappings): {end - start:.2f} s.")


def run():
    parser = argparse.ArgumentParser(
        description="Generate SIFTS mappings between a structure and sequence database."
    )
    parser.add_argument(
        "-i", "--cif-file", required=True,
        help="Path to a PDBx/mmCIF file or a text file listing multiple CIF paths."
    )
    parser.add_argument(
        "-od", "--output-dir", required=True,
        help="Directory where all results will be written."
    )
    parser.add_argument(
        "-db", "--db-file", required=True,
        help="Path to the preformatted sequence database (MMseqs or BLAST)."
    )
    parser.add_argument(
        "-t", "--tool", default="mmseqs",
        help="Alignment tool to use ('mmseqs' or 'blastp')."
    )
    parser.add_argument(
        "-ucsv", "--unp-csv-file", default=None,
        help="Path to the csv file containing accession info: accession, provenance, pdb_xref, annotation_score."
    )
    parser.add_argument(
        "-threads", "--threads", type=int, default=1,
        help="Number of threads to use for parsing and searches."
    )
    parser.add_argument(
        "-bs", "--batch-size", type=int, default=100000,
        help="Number of CIF files to process per batch (default: 100000)."
    )

    args = parser.parse_args()
    logger.info(vars(args))

    pipeline = SiftsGlobalMappings(
        cif_file=args.cif_file,
        out_dir=args.output_dir,
        db_file=args.db_file,
        unp_csv=args.unp_csv_file,
        tool=args.tool,
        threads=args.threads,
        batch_size=args.batch_size,
    )
    pipeline.process()


if __name__ == "__main__":
    run()
