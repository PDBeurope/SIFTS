#!/usr/bin/env python3

import argparse
import subprocess
import shutil
from typing import Union, Optional, List
from pathlib import Path

from pdbe_sifts.base.log import logger
from pdbe_sifts.global_mappings.base_alignment_search import AlignmentSearch

class BlastP(AlignmentSearch):
    def __init__(
        self,
        query_path: Union[str, Path],
        target_path: Union[str, Path],
        output_path: Union[str, Path],
        outfmt: str = "6 qseqid sseqid length mismatch qstart qend sstart send evalue bitscore qseq sseq qlen staxid pident qcovs",
        evalue: float = 10.0,
        threads: int = 1,
        extra_args: Optional[List[str]] = None,
    ):
        super().__init__(query_path, target_path, output_path)
        self.outfmt = outfmt
        self.evalue = evalue
        self.threads = threads
        self.extra_args = extra_args or []

    def _process(self, extra_args: list = None):
        if shutil.which("blastp") is None:
            raise FileNotFoundError("blastp command not found. Please install BLAST+.")

        cmd = [
            "blastp",
            "-query", str(self.query_path),
            "-db", str(self.target_path),
            "-out", str(self.output_path),
            "-outfmt", str(self.outfmt),
            "-evalue", str(self.evalue),
            "-num_threads", str(self.threads),
        ]

        if self.extra_args:
            cmd += self.extra_args

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"blastp execution failed with exit code {e.returncode}")
            raise

def run():
    parser = argparse.ArgumentParser(
        description="Run a blastp search against a blast database."
    )

    parser.add_argument(
        "-query", "--query-path",
        required=True,
        help="Path to the input fasta file.",
    )
    parser.add_argument(
        "-target", "--target-path",
        required=True,
        help="Path to the location of the target blast database.",
    )
    parser.add_argument(
        "-o", "--output-path",
        required=True,
        help="Path to the file where to save the results.",
    )
    parser.add_argument(
        "-outfmt", "--outfmt",
        type=int,
        default=15,
        help="Format of the results. Default: 15 (JSON). See blastp -help for options.",
    )
    parser.add_argument(
        "-eval", "--e-value",
        type=float,
        default=10.0,
        help="Expectation value threshold for saving hits. Default = 10.0",
    )
    parser.add_argument(
        "-threads", "--threads",
        type=int,
        default=1,
        help="Number of threads to use.",
    )
    parser.add_argument(
        "--extra-args",
        nargs=argparse.REMAINDER,
        help="Extra arguments to pass to blastp (e.g., --max_target_seqs 5).",
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
        extra_args=args.extra_args,
    )
    blast_p.run()


if __name__ == "__main__":
    run()