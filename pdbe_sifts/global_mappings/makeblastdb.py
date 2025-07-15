#!/usr/bin/env python3

import argparse
import subprocess
from typing import Union, Optional, List
from pathlib import Path

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.utils import parse_extra_args
from pdbe_sifts.global_mappings.database import ToolDatabase

class MakeBlastDb(ToolDatabase):
    def __init__(
        self,
        fasta_path: Union[str, Path],
        output_path: Union[str, Path],
        tax_id_map: Union[str, Path],
        **kwargs):

        super().__init__(fasta_path, output_path)
        self.kwargs = kwargs

    def _process(self):
        cmd = [
            "makeblastdb",
            "-in", str(self.input_path),
            "-input_type", 'fasta',
            "-dbtype", 'prot',
            "-parse_seqids",
            "-out", str(self.output_path),
            "-taxid_map", str(self.tax_id_map),
        ]

        for key, value in self.kwargs.items():
            flag = f"-{key}"
            if flag not in cmd:
                cmd.append(flag)
                if value is not True or value is not False:
                    cmd.append(str(value))

        subprocess.run(cmd, check=True)

def run():
    parser = argparse.ArgumentParser(
        description="Creation of a blast formated database with makeblastdb command line tool."
    )
    parser.add_argument(
        "-in",
        "--fasta-file",
        required=True,
        help="Base location for the input fasta file. The file should contains at list one sequence.",
    )
    parser.add_argument(
        "-out",
        "--output-db",
        required=True,
        help="Base location where to write the output blast formated database.",
    )
    parser.add_argument(
        "-tax_map",
        "--taxid-map",
        required=True,
        help="Path to the text file containing sequence/tax_id mappings. One row, one mapping separated by tab or space.",
    )
    parser.add_argument(
        "--extra-args",
        nargs=argparse.REMAINDER,
        help="Extra arguments to pass to makeblastdb (e.g., -blastdb_version 4).",
    )

    args = parser.parse_args()
    kwargs = parse_extra_args(args.extra_args)
    logger.info(vars(args))
    blast_db = MakeBlastDb(
        args.fasta_file,
        args.output_db,
        **kwargs,
    )
    blast_db.run()

if __name__ == "__main__":
    run()
