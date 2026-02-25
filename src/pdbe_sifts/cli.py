import argparse
from pathlib import Path

from pdbe_sifts.config import init_config, load_config, _USER_CONFIG_FILE
from pdbe_sifts.sifts_global_mappings import SiftsGlobalMappings

conf = load_config()


def main():
    parser = argparse.ArgumentParser(
        prog="pdbe_sifts",
        description="PDBe SIFTS mapping pipeline"
    )
    subparsers = parser.add_subparsers(dest="command")

    # INIT — copies the YAML template to user config dir
    init_parser = subparsers.add_parser(
        "init",
        help="Copy the config template to ~/.config/pdbe_sifts/config.yaml"
    )
    init_parser.add_argument(
        "-dest", type=Path, default=None,
        help="Custom destination path (default: ~/.config/pdbe_sifts/config.yaml)"
    )
    init_parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing config file."
    )

    # SHOW — display resolved config
    show_parser = subparsers.add_parser(
        "show",
        help="Display the resolved configuration."
    )
    show_parser.add_argument(
        "--config", type=Path, default=None,
        help="Path to a custom config file."
    )

    # RUN sifts_global_mappings
    global_parser = subparsers.add_parser(
        "global_mappings",
        help="Runs alignment and scoring function."
    )
    global_parser.add_argument(
        "-config", "--config", required=False,
        help="Path to a custom config file."
    )
    global_parser.add_argument(
        "-i", "--cif-file", required=True,
        help="Path to a PDBx/mmCIF file or a text file listing multiple CIF paths."
    )
    global_parser.add_argument(
        "-od", "--output-dir", required=False, default=conf.user.base_dir,
        help="Directory where all results will be written."
    )
    global_parser.add_argument(
        "-db", "--db-file", required=False, default=conf.user.target_db,
        help="Path to the preformatted sequence database (MMseqs or BLAST)."
    )
    global_parser.add_argument(
        "-t", "--tool", default="mmseqs",
        help="Alignment tool to use ('mmseqs' or 'blastp')."
    )
    global_parser.add_argument(
        "-ucsv", "--unp-csv-file", default=None,
        help="Path to the csv file containing accession info: accession, provenance, pdb_xref, annotation_score."
    )
    global_parser.add_argument(
        "-threads", "--threads", type=int, default=1,
        help="Number of threads to use for parsing and searches."
    )
    global_parser.add_argument(
        "-bs", "--batch-size", type=int, default=100000,
        help="Number of CIF files to process per batch (default: 100000)."
    )


    args = parser.parse_args()

    if args.command == "init":
        if args.force and _USER_CONFIG_FILE.exists():
            _USER_CONFIG_FILE.unlink()
        init_config(dest=args.dest)

    elif args.command == "show":
        cfg = load_config(args.config)
        print(cfg)
    
    elif args.command == "global_mappings":
        if args.config:
            cfg = load_config(args.config)
        gb_m = SiftsGlobalMappings(
                cif_file=args.cif_file,
                out_dir=args.output_dir,
                db_file=args.db_file,
                unp_csv=args.unp_csv_file,
                tool=args.tool,
                threads=args.threads,
                batch_size=args.batch_size,
            )
        gb_m.process()
        



    # elif args.command == "run":
    #     cfg = load_config(args.config)
    #     if args.overrides:
    #         from omegaconf import OmegaConf
    #         cfg = OmegaConf.merge(cfg, OmegaConf.from_dotlist(args.overrides))
        # ... pipeline logic

    else:
        parser.print_help()
