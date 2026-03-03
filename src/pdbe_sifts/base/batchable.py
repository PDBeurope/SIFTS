"""Batchable — Base class for tasks that process entries one at a time or in parallel.

Usage:

    class MyTask(Batchable):
        workers = 4
        entry_file_path = "/path/to/entries.txt"
        failure_threshold = 0.01

        def worker_setup(self):
            # Called once per worker process before it processes its entries.
            # Put expensive/non-picklable initialization here (e.g. DB connections).
            self.db = open_connection()

        def worker_teardown(self):
            self.db.close()

        def process_entry(self, entry_id: str):
            # Called for every entry in the worker process.
            pass

    # From CLI (in a run() function or __main__):
    custom_args, remaining = parser.parse_known_args()
    task = MyTask(custom_args.foo, custom_args.bar)
    task.main(remaining)   # handles 'single' or 'batch' subcommand

    # Programmatically:
    task = MyTask(...)
    task.run_single("1abc")
    task.run_batch(["1abc", "2xyz", "3pqr"])

Lifecycle — single mode:
    setup() → worker_setup() → process_entry() → worker_teardown() → teardown()

Lifecycle — batch mode:
    setup()
      → [parallel workers: worker_setup() → process_entry()* → worker_teardown()]
      → [retry if enabled: same worker lifecycle for failed entries]
    → teardown()
"""

import argparse
import logging
from abc import ABC, abstractmethod
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from pdbe_sifts.base.exceptions import TooManyFailedEntries
from pdbe_sifts.base.log import logger


# ---------------------------------------------------------------------------
# Module-level worker functions (must be top-level for pickling)
# ---------------------------------------------------------------------------

def _run_chunk(obj: "Batchable", entries: list[str]) -> list[str]:
    """Run a list of entries inside a worker process.

    Calls worker_setup() once, then process_entry() for each entry,
    then worker_teardown() once. Returns list of failed entry IDs.
    """
    failed: list[str] = []
    obj.worker_setup()
    try:
        for entry_id in entries:
            try:
                if obj.entry_timeout_secs and obj.entry_timeout_secs > 0:
                    _run_with_timeout(obj, entry_id, obj.entry_timeout_secs)
                else:
                    obj.process_entry(entry_id)
            except Exception:
                logger.error(f"Failed [{entry_id}]:", exc_info=True)
                failed.append(entry_id)
    finally:
        obj.worker_teardown()
    return failed


def _run_with_timeout(obj: "Batchable", entry_id: str, timeout_secs: int) -> None:
    """Run process_entry in a subprocess with a hard timeout."""
    from multiprocessing import Process

    p = Process(target=obj.process_entry, args=(entry_id,))
    p.start()
    p.join(timeout_secs)
    if p.is_alive():
        p.kill()
        p.join()
        raise TimeoutError(
            f"Entry {entry_id!r} timed out after {timeout_secs}s"
        )


def _split(items: list, n: int) -> list[list]:
    """Split *items* into at most *n* non-empty chunks."""
    if n <= 0 or not items:
        return [list(items)]
    k, rem = divmod(len(items), n)
    chunks: list[list] = []
    i = 0
    for _ in range(n):
        size = k + (1 if rem > 0 else 0)
        rem -= 1
        chunk = items[i : i + size]
        if chunk:
            chunks.append(chunk)
        i += size
    return chunks


def _parse_entry_list(entry_list: str) -> list[str]:
    """Return entry IDs from a file path or a comma-separated string."""
    path = Path(entry_list)
    if path.exists():
        entries = [
            line.strip()
            for line in path.read_text().splitlines()
            if line.strip()
        ]
        logger.info(f"Loaded {len(entries)} entries from {path}")
        return entries
    items = [e.strip() for e in entry_list.split(",") if e.strip()]
    if not items:
        raise ValueError(f"No entries found in: {entry_list!r}")
    return items


# ---------------------------------------------------------------------------
# Batchable base class
# ---------------------------------------------------------------------------

class Batchable(ABC):
    """Base class for tasks that run on individual entries.

    Subclasses **must** implement :meth:`process_entry`.
    All other methods are optional lifecycle hooks.

    Attributes:
        workers: Number of parallel worker processes in batch mode.
        entry_file_path: Default file to load entries from (one ID per line).
            Used when ``--list`` is not supplied on the CLI.
        failure_threshold: Maximum allowed failures before raising.
            A value < 1 is treated as a fraction of total entries (e.g. ``0.01``
            = 1 %); a value >= 1 is an absolute count.
        retry_on_failure: Retry failed entries once after the initial run.
        entry_timeout_secs: Per-entry timeout in seconds. ``None`` = no timeout.

    Note on DuckDB (or any single-writer resource):
        DuckDB does not allow concurrent write connections from multiple
        processes. When ``dbmode=True`` (or equivalent), pass
        ``--workers 1`` to the CLI to avoid locking errors.
    """

    # --- Configuration: override in subclasses ---
    workers: int = 4
    entry_file_path: Optional[str] = None
    failure_threshold: float = 0.0
    retry_on_failure: bool = True
    entry_timeout_secs: Optional[int] = None

    # --- Lifecycle hooks (all optional) ---

    def setup(self) -> None:
        """Called once in the **main process** before any processing starts."""

    def teardown(self) -> None:
        """Called once in the **main process** after all processing is done."""

    def worker_setup(self) -> None:
        """Called once per **worker process** before it processes its entries.

        Use this for expensive or non-picklable initialization such as database
        connections or loading large data structures.
        """

    def worker_teardown(self) -> None:
        """Called once per **worker process** after all its entries are processed.

        Use this to clean up resources allocated in :meth:`worker_setup`.
        """

    # --- Abstract method ---

    @abstractmethod
    def process_entry(self, entry_id: str) -> None:
        """Process a single entry.

        Args:
            entry_id: The entry identifier (e.g. PDB ID, UniProt accession).
        """

    # --- Public run API ---

    def run_single(self, entry_id: str) -> None:
        """Process a single entry with full lifecycle hooks."""
        self.setup()
        try:
            _run_chunk(self, [entry_id])
        finally:
            self.teardown()

    def run_batch(
        self,
        entries: list[str],
        workers: Optional[int] = None,
        retry: Optional[bool] = None,
    ) -> None:
        """Process multiple entries in parallel.

        Args:
            entries: List of entry IDs to process.
            workers: Number of parallel worker processes.
                Defaults to :attr:`workers`.
            retry: Whether to retry failed entries.
                Defaults to :attr:`retry_on_failure`.
        """
        n_workers = workers or self.workers
        should_retry = retry if retry is not None else self.retry_on_failure

        self.setup()
        try:
            failed = self._run_parallel(entries, n_workers)

            if failed and should_retry:
                logger.info(f"Retrying {len(failed)} failed entries...")
                failed = self._run_parallel(failed, min(n_workers, len(failed)))

            if failed:
                self._log_failed(failed)
            self._check_failures(failed, len(entries))
        finally:
            self.teardown()

    def load_entries(self) -> list[str]:
        """Load entries from :attr:`entry_file_path`.

        Override this for custom entry-loading logic.

        Returns:
            List of entry IDs.

        Raises:
            ValueError: If ``entry_file_path`` is not set or the file is empty.
            FileNotFoundError: If the entry file does not exist.
        """
        if not self.entry_file_path:
            raise ValueError(
                "No entry_file_path configured. "
                "Set entry_file_path on the class or pass --list."
            )
        path = Path(self.entry_file_path)
        if not path.exists():
            raise FileNotFoundError(f"Entry file not found: {path}")
        entries = [
            line.strip()
            for line in path.read_text().splitlines()
            if line.strip()
        ]
        if not entries:
            raise ValueError(f"Entry file is empty: {path}")
        logger.info(f"Loaded {len(entries)} entries from {path}")
        return entries

    def add_arguments(self, parser: argparse.ArgumentParser) -> None:
        """Override to add custom CLI arguments.

        This is called with the **root** ArgumentParser before subparsers are
        added, so arguments defined here are available in both ``single`` and
        ``batch`` subcommands.

        Args:
            parser: The root :class:`argparse.ArgumentParser`.
        """

    # --- CLI entry point ---

    def main(self, argv: Optional[list[str]] = None) -> None:
        """Parse CLI arguments and dispatch to :meth:`run_single` or :meth:`run_batch`.

        Subcommands::

            single  --entry <id>
            batch   --list <file_or_csv> [--workers N] [--no-retry]
                    [--failure-threshold F] [--timeout S]

        Args:
            argv: Argument list to parse. Uses ``sys.argv[1:]`` when ``None``.
        """
        root = argparse.ArgumentParser(description=self.__class__.__name__)
        self.add_arguments(root)
        sub = root.add_subparsers(dest="mode", required=True)

        # single subcommand
        single_p = sub.add_parser("single", help="Process one entry.")
        single_p.add_argument(
            "--entry", "-e", required=True, help="Entry ID to process."
        )

        # batch subcommand
        batch_p = sub.add_parser("batch", help="Process multiple entries in parallel.")
        batch_p.add_argument(
            "--list",
            "-l",
            dest="entry_list",
            help=(
                "Path to a file of entry IDs (one per line) or a comma-separated "
                "list. Uses entry_file_path if not provided."
            ),
        )
        batch_p.add_argument(
            "--workers",
            "-w",
            type=int,
            default=self.workers,
            help=f"Number of parallel worker processes (default: {self.workers}).",
        )
        batch_p.add_argument(
            "--no-retry",
            action="store_true",
            help="Disable automatic retry for failed entries.",
        )
        batch_p.add_argument(
            "--failure-threshold",
            type=float,
            default=self.failure_threshold,
            dest="failure_threshold",
            help=(
                "Max allowed failure ratio (<1) or count (>=1). "
                "Default: %(default)s."
            ),
        )
        batch_p.add_argument(
            "--timeout",
            type=int,
            default=self.entry_timeout_secs,
            dest="entry_timeout_secs",
            help="Per-entry timeout in seconds. No timeout if not set.",
        )

        args = root.parse_args(argv)
        logger.info(vars(args))

        if args.mode == "single":
            self.run_single(args.entry)
        else:
            self.failure_threshold = args.failure_threshold
            self.entry_timeout_secs = args.entry_timeout_secs
            entries = (
                _parse_entry_list(args.entry_list)
                if args.entry_list
                else self.load_entries()
            )
            self.run_batch(entries, workers=args.workers, retry=not args.no_retry)

    # --- Internal helpers ---

    def _run_parallel(self, entries: list[str], n_workers: int) -> list[str]:
        """Distribute entries across worker processes; return list of failed IDs."""
        chunks = _split(entries, n_workers)
        failed: list[str] = []

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(_run_chunk, self, chunk): chunk
                for chunk in chunks
            }
            for future in as_completed(futures):
                chunk = futures[future]
                try:
                    chunk_failed = future.result()
                    failed.extend(chunk_failed)
                except Exception:
                    logger.error("Worker process crashed:", exc_info=True)
                    failed.extend(chunk)

        return failed

    def _check_failures(self, failed: list[str], total: int) -> None:
        if not failed:
            logger.info("All entries processed successfully.")
            return
        if self.failure_threshold <= 0:
            return
        allowed = (
            int(self.failure_threshold * total)
            if self.failure_threshold < 1
            else int(self.failure_threshold)
        )
        logger.debug(f"Failures: {len(failed)}, allowed: {allowed}")
        if len(failed) > allowed:
            raise TooManyFailedEntries(
                f"{len(failed)} failures exceed threshold "
                f"({self.failure_threshold} → {allowed} entries)"
            )

    def _log_failed(self, failed: list[str]) -> None:
        count = min(5, len(failed))
        logger.warning(
            f"{len(failed)} entries failed. First {count}: {failed[:count]}"
        )
