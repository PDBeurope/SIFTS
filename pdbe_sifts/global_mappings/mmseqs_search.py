#!/usr/bin/env python3

import argparse

from typing import Union
from pathlib import Path

from pymmseqs.config.easy_search_config import EasySearchConfig
from pdbe_sifts.base.log import logger
from pdbe_sifts.global_mappings.base_alignment_search import AlignmentSearch

class MmSearch(AlignmentSearch):
    def __init__(
        self,
        query_path: Union[str, Path],
        target_path: Union[str, Path],
        output_path: Union[str, Path],
        outtmp_path: Union[str, Path],
        threads: int,
        **kwargs,
    ):
        """Process the search of a query against a sequence database.

        Args:
            query_path (path to file): path to the query fasta file.
            target_path (path to file): path to the DB file.
            output_path (path to file): path to the result file.
            outtmp_path (path to folder): path to the tmp folder needed by mmseqs.
            kwargs: any possible argument allowed by mmseqs easy_search.
        """

        super().__init__(query_path, target_path, output_path)
        self.outtmp_path = outtmp_path
        self.search = None
        self.format_string = 'query,target,alnlen,mismatch,qstart,qend,tstart,tend,evalue,bits,qaln,taln,qlen,taxid,qheader'
        self.threads = threads
        self.easy_search_config_kwargs = kwargs

    def _process(self):
        result = EasySearchConfig(self.query_path,
                                  self.target_path,
                                  self.output_path,
                                  self.outtmp_path,
                                  format_mode=4,
                                  a=True,
                                  alignment_mode = 3,
                                  format_output=self.format_string,
                                  v=3,
                                  threads = self.threads,)
                                #   **self.easy_search_config_kwargs)
        result.run()


def run():
    parser = argparse.ArgumentParser(
        description="Run a mmseqs easy_search against a mmseqs database."
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
        help="Base location of the result file.",
    )

    parser.add_argument(
        "-tmp",
        "--outtmp-path",
        required=True,
        help="Base location of the tmp folder.",
    )
    parser.add_argument(
        "-threads",
        "--threads",
        type=int,
        required=True,
        help="Number of threads to use for the search.",
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