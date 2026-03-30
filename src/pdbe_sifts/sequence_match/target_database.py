#!/usr/bin/env python3

"""
Module for creating a reference target database for sequence search tools.

This module defines the `TargetDb` class, a concrete implementation of `ToolDatabase`
used to build searchable sequence databases for tools such as **MMseqs2** and **BLAST**.
It provides a command-line interface (`run()`) for creating the database directly
from a FASTA file.

Usage (CLI example):
    python target_db.py -i input.fasta -o output_path --tool mmseqs

Classes:
    TargetDb: Concrete subclass of ToolDatabase implementing database creation for MMseqs2 or BLAST.
Functions:
    run(): Command-line entry point for creating a target database.
"""

import argparse
from pathlib import Path

from pymmseqs.commands import createdb
from pymmseqs.config.createindex_config import CreateIndexConfig
from pymmseqs.config.createtaxdb_config import CreateTaxDBConfig

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.paths import get_conf_mmseqs_index_subset
from pdbe_sifts.sequence_match.database import ToolDatabase
from pdbe_sifts.sequence_match.makeblastdb import MakeBlastDb


class TargetDb(ToolDatabase):
    """
    Concrete implementation of ToolDatabase for MMseqs2 or BLAST.

    This class handles the creation of a reference target database that
    alignment tools can search against. It supports both MMseqs2 and
    BLAST, invoking either `mmseqs createdb` and `createindex`, or
    `makeblastdb`, depending on the chosen tool.

    Args:
        input_path (str | Path): Path to the input FASTA file containing sequences.
        output_path (str | Path): Path where output database files will be saved.
        tax_mapping_file (str | Path): File mapping sequence identifiers to taxonomic identifiers.
        tool (str): The tool to use for database creation ('mmseqs' or 'blastp'). Defaults to 'mmseqs'.
        threads (int): Number of CPU threads to use. Defaults to 1.
        **kwargs: Additional keyword arguments passed to the database creation commands.
    """

    def __init__(
        self,
        input_path: str | Path,
        output_path: str | Path,
        tax_mapping_file: str | Path,
        tool: str = "mmseqs",
        threads: int = 1,
        **kwargs,
    ):
        """Initialize the TargetDb instance with input/output paths and configuration."""
        super().__init__(input_path, output_path)
        self.tool = tool
        self.target_db = None
        self.tax_mapping_file = Path(tax_mapping_file).resolve()
        self.threads = threads
        self.db_config_kwargs = kwargs

    def _process(self) -> None:
        """Run the appropriate database creation tool based on the configuration."""
        match self.tool:
            case "mmseqs":
                self.target_db = createdb(
                    self.input_path, self.output_path, **self.db_config_kwargs
                )
                tmp_fold = Path(self.output_path).parent / "tmp"
                tmp_fold.mkdir(parents=True, exist_ok=True)
                CreateTaxDBConfig(
                    sequence_db=self.target_db.to_path(),
                    tmp_dir=tmp_fold,
                    tax_mapping_file=self.tax_mapping_file,
                    threads=self.threads,
                ).run()
                CreateIndexConfig(
                    self.target_db.to_path(),
                    tmp_dir=tmp_fold,
                    threads=self.threads,
                    index_subset=get_conf_mmseqs_index_subset(),
                ).run()
            case "blastp":
                self.target_db = MakeBlastDb(
                    self.input_path, self.output_path, self.tax_mapping_file
                )
                self.target_db._process()


def run() -> None:
    """Command-line entry point for creating a MMseqs2 or BLAST target database."""
    parser = argparse.ArgumentParser(
        description="Create a MMseqs2 or BLAST target database from a FASTA file."
    )

    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        help="Path to the input FASTA file (with at least one sequence).",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=True,
        help="Path where the target database files will be saved.",
    )
    parser.add_argument(
        "-t",
        "--tax-mapping-file",
        required=True,
        help="File mapping sequence identifiers to taxonomic identifiers.",
    )
    parser.add_argument(
        "--tool",
        default="mmseqs",
        help="Tool to use for database creation ('mmseqs' or 'blastp'). Default: mmseqs.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of CPU threads to use. Default: 1.",
    )
    args = parser.parse_args()

    logger.info(vars(args))
    target_db = TargetDb(
        input_path=args.input_file,
        output_path=args.output_path,
        tax_mapping_file=args.tax_mapping_file,
        tool=args.tool,
        threads=args.threads,
    )
    target_db.run()


if __name__ == "__main__":
    run()
