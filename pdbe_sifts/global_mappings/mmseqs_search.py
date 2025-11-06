#!/usr/bin/env python3
"""
Module for running MMseqs2 easy_search.

This module defines the `MmSearch` class, a concrete implementation of `AlignmentSearch`
for performing MMseqs2 sequence-to-database searches. It wraps the `EasySearchConfig`
interface to configure and execute MMseqs2 with controlled parameters.

Classes:
    MmSearch: Concrete subclass of AlignmentSearch implementing MMseqs2 easy_search.
Functions:
    run(): Command-line interface for running an MMseqs2 search.
"""

import argparse
from typing import Union
from pathlib import Path

from pymmseqs.config.easy_search_config import EasySearchConfig
from pdbe_sifts.base.log import logger
from pdbe_sifts.global_mappings.base_alignment_search import AlignmentSearch

class MmSearch(AlignmentSearch):
    """
    Concrete implementation of AlignmentSearch using MMseqs2 easy_search.

    This class runs an MMseqs2 search comparing a query FASTA file against
    a pre-built MMseqs2 target database. The results are saved in a tabular format.

    Args:
        query_path (str | Path): Path to the query FASTA file.
        target_path (str | Path): Path to the MMseqs2 target database.
        output_path (str | Path): Path to the output results file.
        outtmp_path (str | Path): Path to the temporary directory required by MMseqs2.
        threads (int): Number of CPU threads to use for the search.

    Attributes:
        outtmp_path (Path): Temporary folder for intermediate MMseqs2 files.
        search: EasySearchConfig instance used to configure the search.
        format_string (str): Output format for MMseqs2 results.
        threads (int): Number of threads used.
    """
    def __init__(
        self,
        query_path: Union[str, Path],
        target_path: Union[str, Path],
        output_path: Union[str, Path],
        outtmp_path: Union[str, Path],
        threads: int,
        **kwargs,
    ):
        """Initialize the MmSearch instance with query, target, and output paths."""

        super().__init__(query_path, target_path, output_path)
        self.outtmp_path = outtmp_path
        self.search = None
        self.format_string = (
            "query,target,alnlen,mismatch,qstart,qend,tstart,tend,evalue,bits,"
            "qaln,taln,qlen,taxid,qheader,fident,qcov"
        )
        self.threads = threads
        self.easy_search_config_kwargs = kwargs

    def _process(self):
        """Run the MMseqs2 easy_search process."""
        result = EasySearchConfig(self.query_path,
                                  self.target_path,
                                  self.output_path,
                                  self.outtmp_path,
                                  format_mode=4,
                                  a=True,
                                  alignment_mode = 3,
                                  format_output=self.format_string,
                                  v=3,
                                  threads = self.threads,
                                  db_load_mode=2,
                                  s=7.5,
                                  max_seqs = 500,)
                                #   **self.easy_search_config_kwargs)
        result.run()


def run():
    """Command-line interface for running MMseqs2 easy_search."""
    parser = argparse.ArgumentParser(
        description="Run a mmseqs easy_search against a mmseqs database."
    )

    parser.add_argument(
        "-query",
        "--query-path",
        required=True,
        help="Path to the query FASTA file.",
    )
    parser.add_argument(
        "-target",
        "--target-path",
        required=True,
        help="Path to the MMseqs2 target database.",
    )

    parser.add_argument(
        "-o",
        "--output-path",
        required=True,
        help="Path to save the search result file.",
    )

    parser.add_argument(
        "-tmp",
        "--outtmp-path",
        required=True,
        help="Path to the temporary directory required by MMseqs2.",
    )
    parser.add_argument(
        "-threads",
        "--threads",
        type=int,
        required=True,
        help="Number of CPU threads to use for the search.",
    )

    args = parser.parse_args()

    logger.info(vars(args))
    mm_search = MmSearch(
        args.query_path,
        args.target_path,
        args.output_path,
        args.outtmp_path,
        threads=args.threads
    )
    mm_search.run()

if __name__ == "__main__":
    run()