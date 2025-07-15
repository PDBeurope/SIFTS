#!/usr/bin/env python3

import argparse

from pathlib import Path
from typing import Union
from pymmseqs.config.align_config import AlignConfig

from pdbe_sifts.global_mappings.database import ToolDatabase
from pdbe_sifts.base.log import logger
from pdbe_sifts.base.utils import parse_extra_args

class AlignDb(ToolDatabase):
    def __init__(
        self,
        result_db_path: Union[str, Path],
        align_db_path: Union[str, Path],
        query_db_path: Union[str, Path],
        target_db_path: Union[str, Path],
        **kwargs,
    ):
        """
        Generate an alignment database from the result of an MMseqs search using mmseqs align.

        Args:
            result_db_path: Path to the resultDB file from the search (input).
            align_db_path: Path where the output alignDB will be saved.
            query_db_path: Path to the input QueryDB file.
            target_db_path: Path to the input TargetDB file.
            **kwargs: Additional parameters passed to AlignConfig.
        """
        super().__init__(result_db_path, align_db_path)
        self.query_db_path = Path(query_db_path)
        self.target_db_path = Path(target_db_path)
        for path in [self.query_db_path, self.target_db_path]:
            if not path.exists():
                raise FileNotFoundError(f"Required file not found: {path}")
        self.align_config_kwargs = kwargs
        self.align_config_kwargs.setdefault("a", True)
        self.align_config_kwargs.setdefault("alignment_mode", 3)
        self.align_db = None

    def _process(self):
        self.align_db = AlignConfig(
                                    query_db=self.query_db_path,
                                    target_db=self.target_db_path,
                                    result_db=self.input_path,
                                    alignment_db=self.output_path,
                                    **self.align_config_kwargs,
                                )
        self.align_db.run()

def run():
    parser = argparse.ArgumentParser(
        description="Creation of a mmseqs align database from the search result."
    )

    parser.add_argument(
        "-i",
        "--result-db-path",
        required=True,
        help="Base location for the input resultDB file.",
    )
    parser.add_argument(
        "-o",
        "--align-db-path",
        required=True,
        help="Base location where output database files will be saved.",
    )
    parser.add_argument(
        "-q",
        "--query-db-path",
        required=True,
        help="Base location for the input QueryDB file.",
    )
    parser.add_argument(
        "-t",
        "--target-db-path",
        required=True,
        help="Base location for the input TargetDB file.",
    )
    parser.add_argument(
            "--extra-args",
            nargs=argparse.REMAINDER,
            help="Extra arguments passed to mmseqs align (e.g. --extra-args -alignment_mode 3 -a)",
        )

    args = parser.parse_args()
    kwargs = parse_extra_args(args.extra_args or [])
    logger.info(vars(args))
    align_db = AlignDb(
        args.result_db_path,
        args.align_db_path,
        args.query_db_path,
        args.target_db_path,
        **kwargs,
    )
    align_db.run()

if __name__ == "__main__":
    run()
