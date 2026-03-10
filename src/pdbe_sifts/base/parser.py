"""Base parser for all batchable applications.

This module provides a base parser that has common arguments for batchable
applications. It also provides a function to merge the application specific
parser with the base parser and parse the arguments together.

The options provided by the base parser are:
--list: File with list of entries to process or comma-separated list of entries.
--filter: Filter to be used for refinement of input list.
--initial-memory: Initial memory to be used for the application.
--retry-memory: Memory to be used for the application in case of failure.
--no-retry: Do not retry failed entries.
--no-of-jobs: Number of jobs to be run in parallel.
--cores: Number of cores to be used for the application.
--cluster-queue: Name of the cluster queue to be used for the application.
--entry-timeout-secs: Timeout for a single entry in seconds.
--failure-threshold: Number of failures before a job is flagged as failed.
--ignore-batch-failures: Ignore batch failures and exit with success.
--run-local: Run the application locally rather than on the farm.
--main-queue-name: Name of the messaging queue to be used for the application. Other
                    queues (failed, result) will be derived from this name.
                    This name is also used as the job name on the farm
                    and for log directory
--executor: Executor to be used for the application. Default: SIFTS_EXECUTOR env var or conf.executor.type
"""

import argparse
import re
import sys
from pathlib import Path

from pdbe_sifts.base.log import logger
from pdbe_sifts.config import load_config

conf = load_config()
# Defaults
INIT_MEM = 1024
RETRY_MEM = 4096
JOBS = 100
CORES = 1
GPUS = "0"

filter_help = """Filter to be used for refinement of input list.
                 Only entries with an allowed size are going to be
                 processed. The rest is skipped.

                 Usage:
                    '300' -> range(0, 300)
                    '100-300' -> range(100,300)
                    '300+' -> range(300, sys.maxsize)
              """


def filter_type(value: str):
    """Create a filter to be used for data loading.

    Args:
        value (str): Passed application argument

    Raises:
        argparse.ArgumentTypeError: on parsing error.

    Returns:
        range: Entry range to be considered for processing.
    """
    filter_range = None

    if re.fullmatch(r"\d+", value):
        v = int(value)
        if v > 0:
            filter_range = range(0, v)

    elif re.fullmatch(r"\d+\-\d+", value):
        splt = sorted([int(x) for x in value.split("-")])
        filter_range = range(splt[0], splt[1])

    elif re.fullmatch(r"\d+\+", value):
        v = int(value[:-1])
        filter_range = range(v, sys.maxsize)

    if filter_range:
        return filter_range
    else:
        raise argparse.ArgumentTypeError(
            f"Filter argument '{value}' could not be matched."
        )


def parse_with_base_parser(
    app_parser: argparse.ArgumentParser, args=None, namespace=None
) -> argparse.Namespace:
    """Merge the subclass parser into base parsers and parse arguments.

    Args:
        app_parser: Application entry specific parser
        args (List, optional): Use this args rather than `sys.argv`

    Returns:
       Parsed arguments namespace
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="subcommands", dest="subcommand")
    parser_entry = subparsers.add_parser("single", add_help=False, parents=[app_parser])
    parser_entry.add_argument("--entry", required=True, help="Entry ID to process.")

    parser_batch = subparsers.add_parser("batch", add_help=False, parents=[app_parser])
    parser_batch.add_argument(
        "--list",
        help="""
    File with list of entries to process or comma-separated list of entries.
    Must be more more than one item
    """,
    )
    parser_batch.add_argument(
        "--no-of-jobs", type=int, help=f"Number of jobs. Default: {JOBS}"
    )
    parser_batch.add_argument(
        "--cores",
        type=int,
        help=(
            "CPU cores to use. If >1 and run-local is True,"
            f" multiproc is used. Default: {CORES}"
        ),
    )
    parser_batch.add_argument(
        "--gpus",
        help=(
            "Number of gpus to use. If 0, normal node will be used,"
            f" Otherwise, gpu nodes will be used."
            f" It can be passed as gpu:2 or gpu:a100:2"
            f" Default: {GPUS}"
        ),
    )
    parser_batch.add_argument(
        "--initial-memory",
        type=int,
        help=f"Initial memory to request on the farm. Default: {INIT_MEM}",
    )
    parser_batch.add_argument(
        "--retry-memory",
        type=int,
        help=f"Retry memory to request on the farm.Default: {RETRY_MEM}",
    )
    parser_batch.add_argument(
        "--log-dir",
        help="""
        Log directory for LSF jobs. default: conf.location.work.logs_dir/__class__.__name__",
        """,
    )
    parser_batch.add_argument(
        "--run-local",
        action="store_true",
        default=False,
        help="Run locally instead of LSF. Uses multithreading.",
    )
    parser_batch.add_argument(
        "--main-queue-name",
        default=None,
        help="main_queue_name. Autogenerated if not specified",
    )
    parser_batch.add_argument(
        "--pickle-file",
        default=None,
        help="Path to pickle the class to. Shoule be accessible on all workers in farm",
    )
    parser_batch.add_argument(
        "--no-retry",
        default=False,
        action="store_true",
        help="Skip retrying of failed entries",
    )
    parser_batch.add_argument(
        "--cluster-project",
        help="LSF project to use. Ignore if --run-local",
    )
    parser_batch.add_argument(
        "--cluster-queue",
        help="LSF Queue to use. Ignore if --run-local",
    )
    parser_batch.add_argument(
        "--failure-threshold",
        type=float,
        help="Fraction/number of failures allowed before stopping the job. Default: 0.0",
        default=0.0,
    )

    parser_batch.add_argument(
        "--ignore-batch-failures",
        action="store_true",
        help="Ignore failures in batch mode. Overrides --failure-threshold. Default: False",
        default=False,
    )
    parser_batch.add_argument(
        "--executor",
        choices=("lsf", "slurm", "local"),
        help="Executor to use. Default: SIFTS_EXECUTOR env var or conf.executor.type",
    )

    parser_batch.add_argument(
        "--filter", required=False, type=filter_type, help=filter_help
    )
    parser_batch.add_argument(
        "--entry-timeout-secs",
        help="Timeout for an entry in secs, will be marked as failed after this duration",
    )

    args = parser.parse_args(args, namespace)

    if "list" in args and args.list:
        items = args.list.split(",")
        if len(items) == 1 and not Path(args.list).exists():
            parser.error(
                f"""
                File {args.list} does not exist or is not readable.
                If you're trying a single entry, use --single mode
                """
            )
        if len(items) > 1:
            args.list = items

    if not args.subcommand:
        parser.error("subcommand must be specified")
    logger.info(args)
    if args.subcommand == "batch" and args.executor == "local":
        args.run_local = True
    return args
