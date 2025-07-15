#!/usr/bin/env python3

import argparse

from typing import Union
from pathlib import Path
from pymmseqs.commands import createdb
from pdbe_sifts.base.log import logger
from pdbe_sifts.global_mappings.database import ToolDatabase
from pdbe_sifts.base.utils import parse_extra_args

class QueryDb(ToolDatabase):
    def __init__(
        self,
        fasta_path: Union[str, Path],
        output_path: Union[str, Path],
        **kwargs,
    ):
        """Generate the query database to run mmseqs search.

        Args:
            fasta_path (path to file): path to the fasta file (.gz allowed).
            output_path (path to location): location where the output DB will be saved.
        """
        super().__init__(fasta_path, output_path)
        self.query_db = None
        self.query_path = None
        self.query_config_kwargs = kwargs

    def _process(self):
        self.query_db = createdb(self.input_path, self.output_path, **self.query_config_kwargs)
        self.query_path = self.query_db.to_path()

def run():
    parser = argparse.ArgumentParser(
        description="Creation of a mmseqs query database not indexed."
    )

    parser.add_argument(
        "-i",
        "--fasta-path",
        required=True,
        help="Base location for the input fasta file (with at least one sequence).",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=True,
        help="Base location where output database files will be saved.",
    )
    parser.add_argument(
            "--extra-args",
            nargs=argparse.REMAINDER,
            help="Extra arguments passed to mmseqs align (e.g. --extra-args -alignment_mode 3 -a)",
        )

    args = parser.parse_args()
    kwargs = parse_extra_args(args.extra_args or [])
    logger.info(vars(args))
    query_db = QueryDb(
        args.fasta_path,
        args.output_path,
        **kwargs,
    )
    query_db.run()

if __name__ == "__main__":
    run()
