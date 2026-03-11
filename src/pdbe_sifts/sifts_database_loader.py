#!/usr/bin/env python3
"""Load segment/residue CSV files produced by sifts_segments_generation into DuckDB.

Run this once after all batch jobs have completed:

    python sifts_load_segments_db.py \
        -o /path/to/output_dir \
        -db /path/to/sifts_xref.duckdb \
        [--batch-size 1000]
"""

import argparse

import duckdb

from pdbe_sifts.base.log import logger
from pdbe_sifts.database.sifts_db_wrapper import SiftsDB


def run():
    parser = argparse.ArgumentParser(
        description="Bulk-load segment/residue CSVs from sifts_segments_generation into DuckDB."
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Root directory containing *_seg.csv.gz and *_res.csv.gz files.",
    )
    parser.add_argument(
        "-db",
        "--duckdb",
        required=True,
        help="Path to the output DuckDB file (created if it does not exist).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1000,
        help="Number of CSV files per DuckDB transaction (default: 1000).",
    )
    args = parser.parse_args()

    logger.info("Connecting to %s", args.duckdb)
    conn = duckdb.connect(args.duckdb)
    try:
        SiftsDB(conn).bulk_load_from_dir(args.output_dir, args.batch_size)
    finally:
        conn.close()
    logger.info("Done.")


if __name__ == "__main__":
    run()
