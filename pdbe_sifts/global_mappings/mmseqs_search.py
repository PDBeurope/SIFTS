#!/usr/bin/env python3

import argparse

from pymmseqs.commands import search
from pdbe_sifts.base.log import logger
from pdbe_sifts.global_mappings.base_alignment_search import AlignmentSearch

class MmSearch(AlignmentSearch):
    def __init__(
        self,
        query_path,
        target_path,
        output_path,
        outtmp_path,
    ):
        """Process the search of a query against a sequence database.

        Args:
            path_fasta_file (path to file): path to the DB fasta file(s) (.gz allowed).
            out_db (path to location): location where the output DB will be saved.
        """

        super().__init__(query_path, target_path, output_path)
        self.outtmp_path = outtmp_path
        self.search = None

    def _process(self):
        self.search = search(self.query_path, self.target_path, self.output_path, self.outtmp_path)

def run():
    parser = argparse.ArgumentParser(
        description="Run a mmseqs search against a mmseqs database."
    )

    parser.add_argument(
        "-query",
        "--query-path",
        required=True,
        help="Base location for the query file.",
    )
    parser.add_argument(
        "-target",
        "--target-path",
        required=True,
        help="Base location of the targeted database.",
    )

    parser.add_argument(
        "-o",
        "--output-path",
        required=True,
        help="Base location of the result database.",
    )

    parser.add_argument(
        "-tmp",
        "--outtmp-path",
        required=True,
        help="Base location of the tmp folder.",
    )
    args = parser.parse_args()

    logger.info(vars(args))
    mm_search = MmSearch(
        args.query_path,
        args.target_path,
        args.output_path,
        args.outtmp_path,
    )
    mm_search.run()

if __name__ == "__main__":
    run()