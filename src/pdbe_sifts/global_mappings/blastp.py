#!/usr/bin/env python3
"""
Module for running BLASTP searches.

This module defines the `BlastP` class, a concrete implementation of `AlignmentSearch`
that wraps the NCBI BLAST+ `blastp` command-line tool. It provides a consistent
interface and CLI entry point for performing protein sequence alignments against
a pre-built BLAST database.

Classes:
    BlastP: Concrete subclass of AlignmentSearch for running NCBI blastp.
Functions:
    run(): Command-line entry point for performing a BLASTP search.
"""

import argparse
import shutil
import subprocess
from pathlib import Path

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.paths import get_conf_blastp_evalue
from pdbe_sifts.global_mappings.base_alignment_search import AlignmentSearch


class BlastP(AlignmentSearch):
    """
    Concrete implementation of AlignmentSearch using NCBI BLASTP.

    This class runs a `blastp` search comparing a query FASTA file against
    a pre-built BLAST database. It executes the BLAST+ command-line tool
    via `subprocess` and captures logs for reproducibility.

    Args:
        query_path (str | Path): Path to the input FASTA query file.
        target_path (str | Path): Path to the target BLAST database.
        output_path (str | Path): Path to the file where results will be saved.
        outfmt (str, optional): BLAST output format string or integer code. Defaults to a tabular format (6).
        evalue (float, optional): E-value threshold. Defaults to 10.0.
        threads (int, optional): Number of CPU threads to use. Defaults to 1.

    Attributes:
        outfmt (str): Output format setting for BLAST.
        evalue (float): E-value threshold.
        threads (int): Number of threads used for the search.
    """

    def __init__(
        self,
        query_path: str | Path,
        target_path: str | Path,
        output_path: str | Path,
        outfmt: str = "6 qseqid sseqid length mismatch qstart qend sstart send evalue bitscore qseq sseq qlen staxid pident qcovs",
        evalue: float = None,
        threads: int = 1,
    ):
        """Initialize the BlastP instance with search configuration."""
        super().__init__(query_path, target_path, output_path)
        self.outfmt = outfmt
        self.evalue = evalue if evalue is not None else get_conf_blastp_evalue()
        self.threads = threads

    def _process(self):
        """Run the BLASTP search using subprocess."""
        if shutil.which("blastp") is None:
            raise FileNotFoundError("blastp command not found. Please install BLAST+.")

        cmd = [
            "blastp",
            "-query",
            str(self.query_path),
            "-db",
            str(self.target_path),
            "-out",
            str(self.output_path),
            "-outfmt",
            str(self.outfmt),
            "-evalue",
            str(self.evalue),
            "-num_threads",
            str(self.threads),
        ]

        logger.info(f"Running BLASTP: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"blastp execution failed with exit code {e.returncode}")
            raise


def run():
    """Command-line interface for running BLASTP searches."""
    parser = argparse.ArgumentParser(
        description="Run a BLASTP search against a BLAST-formatted database."
    )

    parser.add_argument(
        "-query",
        "--query-path",
        required=True,
        help="Path to the input FASTA file containing protein sequences.",
    )
    parser.add_argument(
        "-target",
        "--target-path",
        required=True,
        help="Path to the BLAST database to search against.",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=True,
        help="Path to save the BLASTP results.",
    )
    parser.add_argument(
        "-outfmt",
        "--outfmt",
        default="6 qseqid sseqid length mismatch qstart qend sstart send evalue bitscore qseq sseq qlen staxid pident qcovs",
        help="BLAST output format (string or integer). Default: tabular format 6.",
    )
    parser.add_argument(
        "-eval",
        "--e-value",
        type=float,
        default=10.0,
        help="E-value threshold for saving hits (default = 10.0).",
    )
    parser.add_argument(
        "-threads",
        "--threads",
        type=int,
        default=1,
        help="Number of CPU threads to use.",
    )

    args = parser.parse_args()
    logger.info(vars(args))

    blast_p = BlastP(
        query_path=args.query_path,
        target_path=args.target_path,
        output_path=args.output_path,
        outfmt=args.outfmt,
        evalue=args.e_value,
        threads=args.threads,
    )
    blast_p.run()


if __name__ == "__main__":
    run()
