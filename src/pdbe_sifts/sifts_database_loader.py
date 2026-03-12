#!/usr/bin/env python3

import argparse

import duckdb

from pdbe_sifts.base.log import logger
from pdbe_sifts.database.sifts_db_wrapper import SiftsDB


def run():
    parser = argparse.ArgumentParser(
        description="Bulk-load segment/residue CSVs from sifts_segments_generation into DuckDB."
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        required=True,
        help="Root directory containing per-entry sifts/ subdirectories with CSV files.",
    )
    parser.add_argument(
        "-d",
        "--duckdb",
        required=True,
        help="Path to the DuckDB file. Entries are read from sifts_xref_segment and reloaded.",
    )
    args = parser.parse_args()

    conn = duckdb.connect(args.duckdb)
    logger.info("Connected to %s", args.duckdb)
    try:
        SiftsDB(conn).bulk_load_from_entries(args.input_dir)
    finally:
        conn.close()
    logger.info("Done.")


if __name__ == "__main__":
    run()
