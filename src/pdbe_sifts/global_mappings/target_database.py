#!/usr/bin/env python3

"""
Module for creating a reference target database for sequence search tools.

This module defines the `TargetDb` class, a concrete implementation of `ToolDatabase`
used to build searchable sequence databases for tools such as **MMseqs2** and **BLAST**.
It provides a command-line interface (`run()`) for creating the database directly
from a FASTA file.

Usage (CLI example):
    python target_db.py -i input.fasta -o output_dir --tool mmseqs

Classes:
    TargetDb: Concrete subclass of ToolDatabase implementing database creation for MMseqs2 or BLAST.
Functions:
    run(): Command-line entry point for creating a target database.
"""

import argparse
from pathlib import Path
from typing import Union

from pymmseqs.commands import createdb, createindex
from pymmseqs.config.createtaxdb_config import CreateTaxDBConfig
from pdbe_sifts.base.log import logger
from pdbe_sifts.global_mappings.makeblastdb import MakeBlastDb
from pdbe_sifts.global_mappings.database import ToolDatabase

class TargetDb(ToolDatabase):
    """
    Concrete implementation of ToolDatabase for MMseqs2 or BLAST.

    This class handles the creation of a reference target database that
    alignment tools can search against. It supports both MMseqs2 and
    BLAST, invoking either `mmseqs createdb` and `createindex`, or
    `makeblastdb`, depending on the chosen tool.

    Args:
        fasta_path (str | Path): Path to the input FASTA file containing sequences.
        output_path (str | Path): Path to the directory where output database files will be saved.
        tax_mapping_file (str | Path): File to map sequence identifier to taxonomical identifier.
        tool (str, optional): The tool to use for database creation ('mmseqs' or 'blastp'). Defaults to 'mmseqs'.
        **kwargs: Additional keyword arguments passed to the database creation commands.
    """
    def __init__(
        self,
        input_path: Union[str, Path],
        output_path: Union[str, Path],
        tax_mapping_file: Union[str, Path],
        tool: Union[str, Path] = 'mmseqs',
        **kwargs,
    ):
        """Initialize the TargetDb instance with input/output paths and configuration."""
        super().__init__(input_path, output_path)
        self.tool = tool
        self.target_db = None
        self.tax_mapping_file = tax_mapping_file
        self.db_config_kwargs = kwargs

    def _process(self):
        """Run the appropriate database creation tool based on the configuration."""
        match self.tool:
            case 'mmseqs':
                self.target_db = createdb(self.input_path, self.output_path,  **self.db_config_kwargs)
                # tmp_fold = Path(self.output_path).parent / 'tmp'
                # tmp_fold.mkdir(parents=True, exist_ok=True)
                # CreateTaxDBConfig(sequence_db=self.target_db.to_path(),
                #                       tmp_dir=tmp_fold,
                #                       tax_mapping_file=self.tax_mapping_file,
                #                       threads=48,).run()
                self.target_db = createindex(self.target_db.to_path())
            case 'blastp':
                self.target_db = MakeBlastDb(self.input_path, self.output_path, self.tax_mapping_file)
                self.target_db._process()

def run():
    """Command-line interface for creating a MMseqs or BLAST target database."""
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
        "-t",
        "--tax-mapping-file",
        required=True,
        help="File to map sequence identifier to taxonomical identifier",
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
        args.tax_mapping_file,
        args.tool,
    )
    target_db.run()

if __name__ == "__main__":
    run()
