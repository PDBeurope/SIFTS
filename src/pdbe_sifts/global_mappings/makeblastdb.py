#!/usr/bin/env python3
"""
Module for creating BLAST-formatted protein sequence databases.

This module defines the `MakeBlastDb` class, a concrete implementation of `ToolDatabase`
that wraps the NCBI BLAST+ `makeblastdb` command-line tool. It provides a consistent
interface and CLI for creating BLAST databases from FASTA files.
"""

import argparse
import subprocess
from typing import Union
from pathlib import Path

from pdbe_sifts.base.log import logger
from pdbe_sifts.global_mappings.database import ToolDatabase


class MakeBlastDb(ToolDatabase):
    """
    Concrete implementation of ToolDatabase using NCBI makeblastdb.

    This class runs `makeblastdb` to create a protein database
    from a FASTA file and an accession-to-taxid mapping file.

    Args:
        input_path (str | Path): Path to the input FASTA file.
        output_path (str | Path): Base path for the output BLAST database.
        tax_id_map (str | Path): Path to the accession-to-taxid mapping file (e.g. "P05067    9606").
    """

    def __init__(
        self,
        input_path: Union[str, Path],
        output_path: Union[str, Path],
        tax_id_map: Union[str, Path],
    ):
        super().__init__(input_path, output_path)
        self.tax_id_map = Path(tax_id_map)

        if not self.input_path.exists():
            raise FileNotFoundError(f"File not found: {self.input_path}")

        if not self.tax_id_map.exists():
            raise FileNotFoundError(f"File not found: {self.tax_id_map}")

    def _process(self):
        """Execute the makeblastdb command to generate the BLAST database."""
        cmd = [
            "makeblastdb",
            "-in", str(self.input_path),
            "-input_type", "fasta",
            "-dbtype", "prot",
            "-parse_seqids",
            "-out", str(self.output_path),
            "-taxid_map", str(self.tax_id_map),
        ]

        logger.info(f"Running makeblastdb: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)


def run():
    """Command-line interface for creating BLAST-formatted databases."""
    parser = argparse.ArgumentParser(
        description="Create a BLAST-formatted protein database using NCBI makeblastdb."
    )
    parser.add_argument(
        "-in", "--fasta-file",
        required=True,
        help="Path to the input FASTA file (must contain at least one protein sequence).",
    )
    parser.add_argument(
        "-out", "--output-db",
        required=True,
        help="Base path for the output BLAST database (no extension).",
    )
    parser.add_argument(
        "-tax_map", "--taxid-map",
        required=True,
        help="Path to the accession-to-taxid mapping file (tab- or space-delimited).",
    )

    args = parser.parse_args()
    logger.info(vars(args))

    blast_db = MakeBlastDb(
        input_path=args.fasta_file,
        output_path=args.output_db,
        tax_id_map=args.taxid_map,
    )

    blast_db.run()


if __name__ == "__main__":
    run()
