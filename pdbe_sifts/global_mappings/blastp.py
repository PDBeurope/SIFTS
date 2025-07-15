#!/usr/bin/env python3

import argparse
import subprocess
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
        outfmt: int = 15,
        evalue: float = 10,
        num_threads: int = 1,
        extra_args: Optional[List[str]] = None,
    ):
        super().__init__(query_path, target_path, output_path)
        self.outfmt = outfmt
        self.evalue = evalue
        self.num_threads = num_threads
        self.extra_args = extra_args or []

    def _process(self, extra_args: list = None):
        cmd = [
            "blastp",
            "-query", str(self.query_path),
            "-db", self.target_path,
            "-out", self.output_path,
            "-outfmt", str(self.outfmt),
            "-evalue", str(self.evalue),
            "-num_threads", str(self.num_threads),
        ] 

        if self.extra_args:
            cmd += self.extra_args

        subprocess.run(cmd, check=True)

def run():
    parser = argparse.ArgumentParser(
        description="Run a blastp search against a blast database."
    )

    parser.add_argument(
        "-query",
        "--query-path",
        required=True,
        help="Path to the input fasta file.",
    )
    parser.add_argument(
        "-target",
        "--target-path",
        required=True,
        help="Path to the location of the target blast database.",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=True,
        help="Path to the file where to save the results.",
    )
    parser.add_argument(
        "-outfmt",
        "--outfmt",
        required=False,
        default='15',
        help="Format of the results. Default: 15 (Json). See blastp -help to know all the options.",
    )
    parser.add_argument(
        "-eval",
        "--e-value",
        required=False,
        default='10',
        help="Expectation value threshold for saving hits. Default = 10",
    )
    parser.add_argument(
        "--extra-args",
        nargs=argparse.REMAINDER,
        help="Extra arguments to pass to blastp (e.g., --max_target_seqs 5).",
    )

    args = parser.parse_args()

    logger.info(vars(args))
    blast_p = BlastP(
        args.query_path,
        args.target_path,
        args.output_path,
        args.outfmt,
        args.e_value,
        extra_args=args.extra_args,
    )
    blast_p.run()

if __name__ == "__main__":
    run()
