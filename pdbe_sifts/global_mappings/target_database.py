#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Union

from pymmseqs.commands import createdb, createindex
from pdbe_sifts.base.log import logger
from pdbe_sifts.global_mappings.makeblastdb import MakeBlastDb
from pdbe_sifts.global_mappings.database import ToolDatabase

class TargetDb(ToolDatabase):
    def __init__(
        self,
        fasta_path: Union[str, Path],
        output_path: Union[str, Path],
        tool: Union[str, Path] = 'mmseqs',
        indexed: bool = False,
        **kwargs,
    ):
        """Generate the reference database against which mmseqs or blast will search. It uses either makeblastdb or mmseqs createdb.

        Args:
            fasta_path (path to file): path to the fasta file containing sequences.
            output_path (path to file): path where the output DB files will be saved.
        """
        super().__init__(fasta_path, output_path)
        self.tool = tool
        self.indexed = True
        self.target_db = None
        self.db_config_kwargs = kwargs

    def _process(self):
        match self.tool:
            case 'mmseqs':
                if self.indexed:
                    self.target_db = createdb(self.input_path, self.output_path,  **self.db_config_kwargs)
                    self.target_db = createindex(self.target_db.to_path())
                else:
                    self.target_db = createdb(self.path_fasta_files, self.out_db, **self.db_config_kwargs)
            case 'blastp':
                self.target_db = MakeBlastDb(self.path_fasta_files, self.out_db)
                self.target_db._process()

def run():
    parser = argparse.ArgumentParser(
        description="Creation of a mmseqs or blast target database."
    )

    parser.add_argument(
        "-i",
        "--fasta-path",
        required=True,
        help="Base location for the fasta file (with at least one sequence).",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=True,
        help="Base location where to saved target database files.",
    )
    parser.add_argument(
        "-tool",
        "--tool",
        required=True,
        help="Tool to use for creating the reference database ('mmseqs' or 'blastp').",
    )
    args = parser.parse_args()

    logger.info(vars(args))
    target_db = TargetDb(
        args.fasta_path,
        args.output_path,
        args.tool,
    )
    target_db._process()

if __name__ == "__main__":
    run()
