#!/usr/bin/env python3
"""
SIFTS search module.

This module performs global sequence–structure mappings for PDB entries
using MMseqs2 or BLASTP. From a fasta file, it runs alignments, and generates a file containing hits (TSV).
"""

import argparse
from timeit import default_timer as timer
from pathlib import Path

from pdbe_sifts.base.log import logger
from pdbe_sifts.global_mappings.mmseqs_search import MmSearch
from pdbe_sifts.global_mappings.blastp import BlastP
from pdbe_sifts.base.utils import get_date, make_path


class SiftsSearchNF:
    """Main class for generating structure–sequence mappings using Nextflow."""

    def __init__(
        self,
        fasta_file: str | Path,
        out_dir: str | Path,
        db_file: str | Path,
        tool: str = "mmseqs",
        threads: int = 1,
    ):
        """
        Initialize SiftsSearchNF (for Nextflow job).

        Args:
            fasta file: Path to a fasta file containing the queries.
            out_dir: Directory where results will be written.
            db_file: Path to the preformatted sequence database (MMseqs or BLAST).
            tool: Alignment tool ('mmseqs' or 'blastp'). Default: 'mmseqs'.
            threads: Number of CPU threads to use.
        """
        self.fasta_file = Path(fasta_file)
        self.out_dir = Path(out_dir)
        self.db_file = Path(db_file)
        self.tool = tool
        self.threads = threads
        self.date = get_date()
        self.result_file_path = {}


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
        logger.info(f"Processing [{self.fasta_file}]")
        start = timer()
        fasta_name = self.fasta_file.name.split('.', 1)[0]
        self.search(fasta_name, self.fasta_file)
        end = timer()


def run():
    parser = argparse.ArgumentParser(
        description="Generate SIFTS mappings between a structure and sequence database."
    )
    parser.add_argument(
        "-i", "--fasta-file", required=True,
        help="Path to a fasta file."
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
        "-threads", "--threads", type=int, default=1,
        help="Number of threads to use for parsing and searches."
    )


    args = parser.parse_args()
    logger.info(vars(args))

    pipeline = SiftsSearchNF(
        fasta_file=args.fasta_file,
        out_dir=args.output_dir,
        db_file=args.db_file,
        tool=args.tool,
        threads=args.threads,
    )
    pipeline.process()


if __name__ == "__main__":
    run()
