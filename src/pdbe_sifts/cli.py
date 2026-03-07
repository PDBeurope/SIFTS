import argparse
from pathlib import Path

from pdbe_sifts.base.log import logger
from pdbe_sifts.config import init_config, load_config, _USER_CONFIG_FILE
from pdbe_sifts.sifts_global_mappings import SiftsGlobalMappings
from pdbe_sifts.global_mappings.target_database import TargetDb

conf = load_config()


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
        "-o", "--output-dir", required=False, default=conf.user.base_dir,
        help="Directory where all results will be written."
    )
    global_parser.add_argument(
        "-d", "--db-file", required=False, default=conf.user.target_db,
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
        help=(
            "Generate SIFTS segments via local sequence alignment (lalign36). "
            "Append 'single --entry <id>' or 'batch --list <file> [--workers N]'."
        ),
    )
    segments_parser.add_argument(
        "-i", "--input-dir",
        default=conf.location.work.data_entry_dir,
        help="Base directory containing mmCIF files.",
    )
    segments_parser.add_argument(
        "-d", "--db-file", required=True,
        help="DuckDB file path.",
    )
    segments_parser.add_argument(
        "-o", "--output-dir",
        default=conf.location.work.data_entry_dir,
        help="Base directory for output CSV files.",
    )
    segments_parser.add_argument(
        "--unp-dir",
        default=conf.cache.uniprot,
        help="Base directory for UniProt cache files.",
    )
    segments_parser.add_argument(
        "--nf90", action="store_true", default=False,
        help="Enable UniRef90 mode (default: False).",
    )
    segments_parser.add_argument(
        "-w", "--write-to-db", action="store_true", default=False,
        help="Also write results to the DuckDB file (default: False).",
    )
    segments_parser.add_argument(
        "-m", "--mapping",
        help="User-defined UniProt accessions, e.g. A:P00963,B:P00963.",
    )
    segments_parser.add_argument(
        "run_args", nargs=argparse.REMAINDER,
        help="Batchable subcommand: 'single --entry <id>' or 'batch --list <file> [--workers N]'.",
    )

    ######### SIFTS → mmCIF
    sifts2mmcif_parser = subparsers.add_parser(
        "sifts2mmcif",
        help="Inject SIFTS mappings into mmCIF files.",
    )
    sifts2mmcif_parser.add_argument(
        "-o", "--output-dir", required=True,
        help="Directory where output mmCIF files will be stored.",
    )
    sifts2mmcif_parser.add_argument(
        "-i", "--input-dir",
        default=conf.location.work.data_entry_dir,
        help="Base directory containing input mmCIF files.",
    )
    sifts2mmcif_parser.add_argument(
        "-s", "--sifts-csv-dir", required=False,
        help="SIFTS CSV base directory. If not given, data is read from the database.",
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
    sifts2mmcif_parser.add_argument(
        "-l", "--delta-sifts-file",
        default=conf.lists.sifts_mapping_changes,
        help="Output file for delta mapping changes (batch mode only).",
    )
    sifts2mmcif_parser.add_argument(
        "--log-dir", required=False, default=None,
        help="Log directory for generating the delta mapping list in teardown().",
    )
    sifts2mmcif_parser.add_argument(
        "run_args", nargs=argparse.REMAINDER,
        help="Batchable subcommand: 'single --entry <id>' or 'batch --list <file> [--workers N]'.",
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
        db_b._process()

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
        sifts_align = SiftsAlign(
            cif_dir=args.input_dir,
            out_dir=args.output_dir,
            file_duckdb=args.db_file,
            unp_dir=args.unp_dir,
            nf90_mode=args.nf90,
            unp_mode=args.mapping,
            dbmode=args.write_to_db,
        )
        sifts_align.main(args.run_args)

    elif args.command == "sifts2mmcif":
        from pdbe_sifts.sifts_to_mmcif.main import ExportSIFTSTommCIF
        obj = ExportSIFTSTommCIF(
            output_path=args.output_dir,
            cif_dir=args.input_dir,
            sifts_csv_dir=args.sifts_csv_dir,
            duckdb_file=args.db_file,
            track_changes=not args.no_track_changes,
            prev_run_dir=args.prev_run_dir,
            delta_file=args.delta_sifts_file,
            log_dir=args.log_dir,
        )
        obj.main(args.run_args)

    else:
        parser.print_help()
