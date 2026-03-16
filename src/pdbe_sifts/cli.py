import argparse
from pathlib import Path

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.paths import (
    get_conf_user_base_dir,
    get_conf_user_target_db,
)
from pdbe_sifts.config import init_config, load_config, _USER_CONFIG_FILE
from pdbe_sifts.sifts_global_mappings import SiftsGlobalMappings
from pdbe_sifts.global_mappings.target_database import TargetDb


def main():
    parser = argparse.ArgumentParser(
        prog="pdbe_sifts",
        description="PDBe SIFTS mapping pipeline"
    )

    # Global flag: applies to every subcommand.
    # Falls back to the SIFTS_LOG_LEVEL environment variable (default: INFO).
    parser.add_argument(
        "--log-level",
        default=None,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        metavar="LEVEL",
        help=(
            "Set the logging verbosity level "
            "(DEBUG, INFO, WARNING, ERROR, CRITICAL). "
            "Overrides the SIFTS_LOG_LEVEL environment variable. "
            "Default: INFO."
        ),
    )

    subparsers = parser.add_subparsers(dest="command")

    #########  INIT — copies the YAML template to user config dir
    init_parser = subparsers.add_parser(
        "init",
        help="Copy the config template to ~/.config/pdbe_sifts/config.yaml"
    )
    init_parser.add_argument(
        "--dest", type=Path, default=None,
        help="Custom destination path (default: ~/.config/pdbe_sifts/config.yaml)"
    )
    init_parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing config file."
    )

    #########  SHOW — display resolved config
    show_parser = subparsers.add_parser(
        "show",
        help="Display the resolved configuration."
    )
    show_parser.add_argument(
        "--config", type=Path, default=None,
        help="Path to a custom config file."
    )

    ######### BUILD DATABASE
    targetbuild_parser = subparsers.add_parser(
        "build_db",
        help="Create a reference target database that alignment tools can search against."
    )
    targetbuild_parser.add_argument(
        "-i", "--input-file", required=True,
        help="Path to the input FASTA file (with at least one sequence)."
    )
    targetbuild_parser.add_argument(
        "-o", "--output-path", required=True,
        help="Path where the target database files will be saved."
    )
    targetbuild_parser.add_argument(
        "-t", "--tax-mapping-file", required=True,
        help="File mapping sequence identifiers to taxonomic identifiers."
    )
    targetbuild_parser.add_argument(
        "--tool", default="mmseqs",
        help="Tool to use for database creation ('mmseqs' or 'blastp'). Default: mmseqs."
    )
    targetbuild_parser.add_argument(
        "--threads", type=int, default=1,
        help="Number of CPU threads to use."
    )

    ######### RUN sifts_global_mappings
    global_parser = subparsers.add_parser(
        "global_mappings",
        help="Run alignment and scoring to generate global SIFTS mappings."
    )
    global_parser.add_argument(
        "-i", "--input-file", required=True,
        help=(
            "Input file: a single mmCIF (.cif / .cif.gz), "
            "a FASTA (.fasta / .fa / .faa), "
            "or a text file listing mmCIF paths (.txt)."
        ),
    )
    global_parser.add_argument(
        "-o", "--output-dir", required=False, default=get_conf_user_base_dir(),
        help="Directory where all results will be written."
    )
    global_parser.add_argument(
        "-d", "--db-file", required=False, default=get_conf_user_target_db(),
        help="Path to the preformatted sequence database (MMseqs or BLAST)."
    )
    global_parser.add_argument(
        "--tool", default="mmseqs",
        help="Alignment tool to use ('mmseqs' or 'blastp'). Default: mmseqs."
    )
    global_parser.add_argument(
        "--unp-csv-file", default=None,
        help="Path to CSV with accession metadata: accession, provenance, pdb_xref, annotation_score."
    )
    global_parser.add_argument(
        "--threads", type=int, default=1,
        help="Number of CPU threads to use for parsing and searches."
    )
    global_parser.add_argument(
        "--batch-size", type=int, default=100000,
        help="Number of CIF files to process per batch when using a .txt list (default: 100000)."
    )

    ######### BUILD FASTA
    fasta_parser = subparsers.add_parser(
        "fasta_build",
        help="Build a query FASTA from mmCIF files."
    )
    fasta_parser.add_argument(
        "-i", "--input-file", required=True,
        help="Input file (.cif, .cif.gz, .fasta, .fa, .faa, or .txt list of CIF paths).",
    )
    fasta_parser.add_argument(
        "-o", "--output-dir", required=True,
        help="Directory where the generated FASTA will be written.",
    )
    fasta_parser.add_argument(
        "--threads", type=int, default=1,
        help="Number of parallel workers for .txt list processing (default: 1).",
    )
    fasta_parser.add_argument(
        "--batch-size", type=int, default=100000,
        help="Number of CIF files per batch for .txt list processing (default: 100000).",
    )

    ######### SEGMENTS GENERATION
    segments_parser = subparsers.add_parser(
        "segments",
        help="Generate SIFTS segments for a single mmCIF entry.",
    )
    segments_parser.add_argument(
        "-i", "--input-cif", required=True,
        help="Input CIF file (.cif / .cif.gz).",
    )
    segments_parser.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory for CSV files.",
    )
    segments_parser.add_argument(
        "-d", "--db-file", required=False, default=None,
        help="DuckDB file path. Optional when -m/--mapping is provided.",
    )
    segments_parser.add_argument(
        "--nf90", action="store_true", default=False,
        help="Enable UniRef90 mode (default: False).",
    )
    segments_parser.add_argument(
        "--no-connectivity", dest="connectivity", action="store_false", default=True,
        help="Disable connectivity mode (default: enabled).",
    )
    segments_parser.add_argument(
        "-m", "--mapping",
        help=(
            "User-defined mapping: UniProt accessions 'A:P00963,B:P00963' "
            "or path to a FASTA file with headers >{auth_asym_id}|{sequence_id}."
        ),
    )
    segments_parser.add_argument(
        "--entry", required=False, default=None,
        help="PDB entry ID. If omitted, derived from _entry.id in the CIF.",
    )

    ######### DB LOAD
    db_load_parser = subparsers.add_parser(
        "db_load",
        help="Bulk-load segment/residue CSVs from segments generation into DuckDB.",
    )
    db_load_parser.add_argument(
        "-i", "--input-dir", required=True,
        help="Root directory containing per-entry sifts/ subdirectories with CSV files.",
    )
    db_load_parser.add_argument(
        "-d", "--duckdb", required=True,
        help="Path to the DuckDB file.",
    )

    ######### UPDATE NCBI
    subparsers.add_parser(
        "update_ncbi",
        help="Update the local NCBI taxonomy database (ete4).",
    )

    ######### SIFTS → mmCIF
    sifts2mmcif_parser = subparsers.add_parser(
        "sifts2mmcif",
        help="Inject SIFTS mappings into a mmCIF file.",
    )
    sifts2mmcif_parser.add_argument(
        "--entry", required=False, default=None,
        help="PDB entry ID. If omitted, derived from _entry.id in the CIF.",
    )
    sifts2mmcif_parser.add_argument(
        "-i", "--input-cif", required=True,
        help="Input CIF file (*_updated.cif.gz).",
    )
    sifts2mmcif_parser.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory where SIFTS mmCIF files will be written.",
    )
    sifts2mmcif_parser.add_argument(
        "-s", "--sifts-csv-dir", required=False,
        help="Flat directory with {entry}_seg.csv.gz / _res.csv.gz. If not given, read from DB.",
    )
    sifts2mmcif_parser.add_argument(
        "-d", "--db-file", required=True,
        help="DuckDB file path.",
    )
    sifts2mmcif_parser.add_argument(
        "-T", "--no-track-changes", action="store_true", default=False,
        help="Disable comparison with previous run to track mapping changes.",
    )
    sifts2mmcif_parser.add_argument(
        "-p", "--prev-run-dir", required=False,
        help="Compare sifts_only.mmcif from this directory for delta tracking.",
    )

    args = parser.parse_args()

    # Apply log level from CLI flag (overrides the env-var default set at import time).
    if args.log_level:
        logger.setLevel(args.log_level)
        try:
            import coloredlogs
            coloredlogs.install(level=args.log_level, logger=logger)
        except ImportError:
            pass

    if args.command == "init":
        if args.force and _USER_CONFIG_FILE.exists():
            _USER_CONFIG_FILE.unlink()
        init_config(dest=args.dest)
        logger.info("Initializing NCBI taxonomy database (first run may download ~70 MB)...")
        try:
            from ete4 import NCBITaxa
            NCBITaxa()
            logger.info("NCBI taxonomy database ready.")
        except Exception as e:
            logger.warning(f"NCBI taxonomy initialization failed: {e}")

    elif args.command == "show":
        cfg = load_config(args.config)
        print(cfg)

    elif args.command == "global_mappings":
        gb_m = SiftsGlobalMappings(
            input_file=args.input_file,
            out_dir=args.output_dir,
            db_file=args.db_file,
            unp_csv=args.unp_csv_file,
            tool=args.tool,
            threads=args.threads,
            batch_size=args.batch_size,
        )
        gb_m.process()

    elif args.command == "build_db":
        db_b = TargetDb(
            input_path=args.input_file,
            output_path=args.output_path,
            tax_mapping_file=args.tax_mapping_file,
            tool=args.tool,
            threads=args.threads,
        )
        db_b.run()

    elif args.command == "fasta_build":
        from pdbe_sifts.sifts_fasta_builder import FastaBuilder
        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = FastaBuilder(
            input_path=args.input_file,
            out_dir=out_dir,
            threads=args.threads,
            batch_size=args.batch_size,
        ).build()
        print(f"FASTA written to: {fasta_path}")

    elif args.command == "segments":
        from pdbe_sifts.sifts_segments_generation import SiftsAlign
        from gemmi import cif as gcif

        if not args.db_file and not args.mapping:
            segments_parser.error("At least one of -d/--db-file or -m/--mapping is required.")

        entry_id = args.entry
        if not entry_id:
            block = gcif.read(str(args.input_cif)).sole_block()
            entry_id = block.find_value("_entry.id").strip('"').lower()

        sifts_align = SiftsAlign(
            args.input_cif,
            args.output_dir,
            args.db_file,
            nf90_mode=args.nf90,
            unp_mode=args.mapping,
            connectivity_mode=args.connectivity,
        )
        sifts_align.process_entry(entry_id)
        if sifts_align.conn:
            sifts_align.conn.close()

    elif args.command == "db_load":
        import duckdb as _duckdb
        from pdbe_sifts.database.sifts_db_wrapper import SiftsDB
        conn = _duckdb.connect(args.duckdb)
        try:
            SiftsDB(conn).bulk_load_from_entries(args.input_dir)
        finally:
            conn.close()
        logger.info("Done.")

    elif args.command == "update_ncbi":
        logger.info("Updating NCBI taxonomy database...")
        from ete4 import NCBITaxa
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()
        logger.info("NCBI taxonomy database updated.")

    elif args.command == "sifts2mmcif":
        from pdbe_sifts.sifts_to_mmcif.main import ExportSIFTSTommCIF
        from gemmi import cif as gcif
        entry_id = args.entry
        if not entry_id:
            block = gcif.read(str(args.input_cif)).sole_block()
            entry_id = block.find_value("_entry.id").strip('"').lower()
        obj = ExportSIFTSTommCIF(
            args.input_cif,
            args.output_dir,
            args.sifts_csv_dir,
            args.db_file,
            not args.no_track_changes,
            args.prev_run_dir,
        )
        try:
            obj.process_entry(entry_id)
        finally:
            obj.conn.close()

    else:
        parser.print_help()
