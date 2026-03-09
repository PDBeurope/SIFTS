"""Batchable — Base class for tasks that process entries one at a time or in parallel.

Supports two execution backends:

  * **Local** (default) — multiprocessing.Pool with imap_unordered, one task
    per entry.  Workers are recycled after ``SIFTS_MAXTASKS_PER_CHILD`` entries
    (default 500) to release accumulated Python heap back to the OS.

  * **Cluster** (LSF / Slurm) — entries are pushed into a shared queue; *N*
    independent cluster jobs each pop-and-process entries until the queue is
    empty.  Requires ``executor.type``, ``queue.*``, and
    ``location.work.pickle_dir`` to be set in config.yaml.

Usage (local mode — unchanged from previous API)::

    class MyTask(Batchable):
        workers = 4
        failure_threshold = 0.01

        def worker_setup(self):
            self.db = open_connection()      # called once per worker process

        def worker_teardown(self):
            self.db.close()

        def process_entry(self, entry_id: str):
            pass

    task = MyTask(...)
    task.main(remaining_argv)   # single / batch subcommand

Usage (cluster mode)::

    pdbe_sifts segments -d hits.duckdb -o ./out batch \\
      --list entries.txt --jobs 50 --memory 8192 --cluster lsf

Lifecycle — single mode:
    setup() → worker_setup() → process_entry() → worker_teardown() → teardown()

Lifecycle — local batch mode:
    setup()
      → [workers: worker_setup() once (via _init_worker), process_entry()* per entry]
    → teardown()

Lifecycle — cluster batch mode:
    setup()
      → push entries to queue, pickle self, submit N cluster jobs
      → each job: before_job_start() → process_entry()* → after_job_end()
        (before_job_start / after_job_end delegate to worker_setup / worker_teardown
         unless overridden)
    → teardown()
"""

import argparse
import logging
import os
import pickle  # nosec
import subprocess
import sys
import traceback
from abc import ABC, abstractmethod
from collections.abc import Iterable
from datetime import datetime
from multiprocessing import Pool, Process
from pathlib import Path
from typing import Optional

from pdbe_sifts.base.exceptions import BatchRunException, TooManyFailedEntries
from pdbe_sifts.base.log import logger


# ---------------------------------------------------------------------------
# Module-level worker state (one copy per worker process, set by _init_worker)
# ---------------------------------------------------------------------------

_WORKER_STATE: Optional["Batchable"] = None


def _init_worker(obj: "Batchable") -> None:
    """Pool initialiser: called once when a worker process is forked.

    Sets the module-level _WORKER_STATE and runs worker_setup().
    Registers worker_teardown() as an atexit handler so it runs on clean exit
    (including when the worker is recycled after maxtasksperchild tasks).
    On SIGKILL (OOM), atexit does not fire — that is acceptable since
    read-only DuckDB connections are safe to abandon.
    """
    import atexit
    global _WORKER_STATE
    _WORKER_STATE = obj
    obj.worker_setup()
    atexit.register(obj.worker_teardown)


def _run_entry(entry_id: str) -> tuple[str, Optional[str]]:
    """Process one entry inside a worker process.

    Returns ``(entry_id, None)`` on success, or ``(entry_id, traceback_str)``
    on failure.  Errors are captured here so they can be reported in the main
    process without raising — this keeps the Pool's result iterator alive for
    subsequent entries.
    """
    try:
        if _WORKER_STATE.entry_timeout_secs and _WORKER_STATE.entry_timeout_secs > 0:
            _run_with_timeout(_WORKER_STATE, entry_id, _WORKER_STATE.entry_timeout_secs)
        else:
            _WORKER_STATE.process_entry(entry_id)
        return (entry_id, None)
    except Exception:
        return (entry_id, traceback.format_exc())


def _run_with_timeout(obj: "Batchable", entry_id: str, timeout_secs: int) -> None:
    """Run process_entry in a subprocess with a hard timeout."""
    p = Process(target=obj.process_entry, args=(entry_id,))
    p.start()
    p.join(timeout_secs)
    if p.is_alive():
        p.kill()
        p.join()
        raise TimeoutError(
            f"Entry {entry_id!r} timed out after {timeout_secs}s"
        )


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

    Local-mode attributes:
        workers: Number of parallel worker processes in local batch mode.
        entry_file_path: Default file to load entries from (one ID per line).
        failure_threshold: Max allowed failures before raising.
            ``< 1`` = fraction of total entries; ``>= 1`` = absolute count.
        retry_on_failure: Retry failed entries once after the initial run.
        entry_timeout_secs: Per-entry timeout in seconds. ``None`` = no timeout.

    Cluster-mode attributes:
        initial_memory: Memory (MB) requested for each initial cluster job.
        retry_memory: Memory (MB) requested when retrying failed entries.
        jobs: Maximum number of parallel cluster jobs.
        cores: CPU cores per cluster job.
        gpus: GPU count (as string) per cluster job.
        cluster_queue: Cluster queue/partition name.
            Falls back to ``conf.farm.cluster_queue`` from config.
        log_dir: Directory for cluster job stdout/stderr logs.
            Falls back to ``conf.executor.log_base_dir`` from config.
        known_exceptions: Set of entry IDs whose failures should be silently
            ignored (logged as warnings only).
        exceptions_file: Path to a file of known-exception entry IDs.
        ignore_batch_failures: If ``True``, do not raise even when
            ``failure_threshold`` is exceeded.
        used_cif_categories: mmCIF categories monitored by this task.
            Used by :meth:`no_used_cif_category_modified`.

    Environment variables:
        SIFTS_MAXTASKS_PER_CHILD: Worker recycling interval in local mode
            (default 500).  Set to ``0`` to disable recycling.
        SIFTS_EXECUTOR: Override ``executor.type`` from config.
        SIFTS_EXECUTOR_LOG_DIR: Override ``executor.log_base_dir`` from config.
        SIFTS_PICKLE_DIR: Override ``location.work.pickle_dir`` from config.
    """

    # --- Local parallel config ---
    workers: int = 4
    entry_file_path: Optional[str] = None
    failure_threshold: float = 0.0
    retry_on_failure: bool = True
    entry_timeout_secs: Optional[int] = None
    failed_entries_file: Optional[Path] = None
    """If set, failed entry IDs are written here (one per line) after each local run."""

    # --- Cluster config ---
    initial_memory: int = 1024
    """Initial memory in MB to request for each cluster job."""

    retry_memory: int = 4096
    """Memory in MB to request when retrying failed entries on the cluster."""

    jobs: int = 100
    """Maximum number of parallel cluster jobs to submit."""

    cores: int = 1
    """CPU cores per cluster job."""

    gpus: str = "0"
    """GPU count per cluster job (as a string)."""

    cluster_queue: str = ""
    """Cluster queue/partition name.  Falls back to conf.farm.cluster_queue."""

    log_dir: str = ""
    """Directory for cluster job stdout/stderr.  Falls back to conf.executor.log_base_dir."""

    known_exceptions: set = set()
    """Entry IDs whose failures are silently ignored (logged as warnings)."""

    exceptions_file: Optional[str] = None
    """Path to a file containing known-exception entry IDs (one per line)."""

    ignore_batch_failures: bool = False
    """If True, never raise TooManyFailedEntries."""

    used_cif_categories: Iterable[str] = set()
    """mmCIF categories monitored by this task (for :meth:`no_used_cif_category_modified`)."""

    # --- Runtime state (set during run) ---
    entries: list[str] = []
    """Current batch of entries (populated in cluster mode)."""

    failed_entries: set = set()
    """Failed entry IDs collected after a cluster batch run."""

    # Internal
    _task_name: str = ""
    _start_timestamp: str = ""
    main_queue: str = ""
    failed_queue: str = ""
    pickle_file: Optional[Path] = None

    # ---------------------------------------------------------------------------
    # Lifecycle hooks (all optional)
    # ---------------------------------------------------------------------------

    def setup(self) -> None:
        """Called once in the **main process** before any processing starts."""

    def teardown(self) -> None:
        """Called once in the **main process** after all processing is done."""

    def worker_setup(self) -> None:
        """Called once per **worker process** before it processes its entries.

        In local mode: called via the Pool initialiser, once per worker fork.
        In cluster mode: called via :meth:`before_job_start` on each compute node.

        Use this for expensive or non-picklable initialisation such as database
        connections or loading large data structures.  Attributes set here are
        **not** pickled back to the main process.
        """

    def worker_teardown(self) -> None:
        """Called once per **worker process** when it exits cleanly.

        In local mode: registered as an atexit handler in each worker.
        In cluster mode: called via :meth:`after_job_end` on each compute node.
        """

    def before_job_start(self) -> None:
        """Called once at the start of each cluster job (ORC compatibility hook).

        Defaults to calling :meth:`worker_setup`.  Override this if you need
        cluster-specific setup that differs from local-mode setup.
        """
        self.worker_setup()

    def after_job_end(self) -> None:
        """Called once at the end of each cluster job (ORC compatibility hook).

        Defaults to calling :meth:`worker_teardown`.
        """
        self.worker_teardown()

    def before_retry(self) -> bool:
        """Called before retrying failed entries on the cluster.

        Override to add pre-retry logic or to conditionally cancel the retry.

        Returns:
            ``True`` to proceed with the retry; ``False`` to skip it.
        """
        return True

    def job_cleanup(self) -> None:
        """Override to perform per-job cleanup on cluster nodes."""

    # ---------------------------------------------------------------------------
    # Abstract method
    # ---------------------------------------------------------------------------

    @abstractmethod
    def process_entry(self, entry_id: str) -> None:
        """Process a single entry.

        Args:
            entry_id: The entry identifier (e.g. PDB ID, UniProt accession).
        """

    # ---------------------------------------------------------------------------
    # Public run API (local mode — unchanged from previous version)
    # ---------------------------------------------------------------------------

    def run_single(self, entry_id: str) -> None:
        """Process a single entry with full lifecycle hooks."""
        self.setup()
        try:
            self.worker_setup()
            try:
                if self.entry_timeout_secs and self.entry_timeout_secs > 0:
                    _run_with_timeout(self, entry_id, self.entry_timeout_secs)
                else:
                    self.process_entry(entry_id)
            except Exception:
                logger.error(f"Failed [{entry_id}]:", exc_info=True)
            finally:
                self.worker_teardown()
        finally:
            self.teardown()

    def run_batch(
        self,
        entries: list[str],
        workers: Optional[int] = None,
        retry: Optional[bool] = None,
    ) -> None:
        """Process multiple entries in parallel using the local multiprocessing Pool.

        This method is the **local** execution path (unchanged from previous
        API).  For cluster submission, use :meth:`run_cluster_batch` or pass
        ``--jobs`` / ``--cluster`` to the CLI which will call it automatically.

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

    # ---------------------------------------------------------------------------
    # Cluster run API
    # ---------------------------------------------------------------------------

    def run_cluster_batch(
        self,
        entries: list[str],
        jobs: Optional[int] = None,
        memory: Optional[int] = None,
        retry_memory: Optional[int] = None,
        cluster: Optional[str] = None,
        cluster_queue: Optional[str] = None,
        queue_name: Optional[str] = None,
        log_dir: Optional[str] = None,
        no_retry: bool = False,
    ) -> None:
        """Submit entries to the cluster via LSF or Slurm.

        Entries are pushed to a shared queue; *jobs* independent cluster jobs
        each pop-and-process entries until the queue is empty.  Failed entries
        are collected in a separate queue and optionally retried with
        higher memory.

        Args:
            entries: List of entry IDs to process.
            jobs: Maximum number of parallel cluster jobs.
                Defaults to :attr:`jobs`.
            memory: Initial memory per job in MB.
                Defaults to :attr:`initial_memory`.
            retry_memory: Memory for retry jobs in MB.
                Defaults to :attr:`retry_memory`.
            cluster: Executor type — ``"lsf"`` or ``"slurm"``.
                Falls back to config / env.
            cluster_queue: Cluster partition/queue name.
                Falls back to :attr:`cluster_queue` or config.
            queue_name: Name used for the work and failed queues.
                Auto-generated from class name + timestamp if not given.
            log_dir: Directory for cluster job logs.
                Falls back to :attr:`log_dir` or config.
            no_retry: If ``True``, skip the automatic retry pass.
        """
        from pdbe_sifts.base.executors import get_executor

        self._task_name = self.__class__.__name__
        self._start_timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

        # Resolve runtime parameters
        n_jobs = jobs or self.jobs
        init_mem = memory or self.initial_memory
        retry_mem = retry_memory or self.retry_memory
        retry_mem = max(init_mem, retry_mem)

        if cluster_queue:
            self.cluster_queue = cluster_queue

        self.main_queue = queue_name or "_".join(
            [self._task_name.lower(), self._start_timestamp]
        )
        self.failed_queue = f"{self.main_queue}-failed"

        resolved_log_dir = (
            log_dir
            or os.getenv("SIFTS_EXECUTOR_LOG_DIR")
            or self.log_dir
            or self._get_config_log_base_dir()
        )
        self.log_dir = str(
            Path(resolved_log_dir, self._task_name.lower(), self._start_timestamp)
        )
        os.makedirs(self.log_dir, exist_ok=True)

        self.executor = get_executor(cluster)
        self.entries = list(entries)

        self.setup()
        try:
            self.failed_entries = set(
                self._submit_batch_internal(n_jobs, init_mem)
            )

            if not no_retry and self.failed_entries and self.before_retry():
                logger.info(f"Retrying {len(self.failed_entries)} failed entries")
                queue = self._get_queue()
                queue.rename(self.failed_queue, self.main_queue)
                retry_jobs = min(len(self.failed_entries), n_jobs)
                self.failed_entries = set(
                    self._submit_batch_internal(retry_jobs, retry_mem)
                )

            if self.failed_entries:
                self.write_failed_entries_to_file()
                self.log_failed_entries()
            else:
                logger.info("All entries processed successfully")

            self._check_failures(list(self.failed_entries), len(entries))

        finally:
            # Remove the pickle file created by _pickle_obj
            if self.pickle_file and Path(self.pickle_file).exists():
                Path(self.pickle_file).unlink()
            self.teardown()
            self._send_summary()

    def _get_config_log_base_dir(self) -> str:
        try:
            from pdbe_sifts.config import load_config
            return load_config().executor.log_base_dir or str(Path.home())
        except Exception:
            return str(Path.home())

    def _get_pickle_dir(self) -> Path:
        env_val = os.getenv("SIFTS_PICKLE_DIR")
        if env_val:
            return Path(env_val)
        try:
            from pdbe_sifts.config import load_config
            cfg = load_config()
            if cfg.location.work.pickle_dir:
                return Path(cfg.location.work.pickle_dir)
        except Exception:
            pass
        return Path.home()

    def _get_queue(self):
        from pdbe_sifts.base.queues.batchable_queue import BatchableQueue
        return BatchableQueue().get_queue()

    def _submit_batch_internal(self, n_jobs: int, memory: int) -> list[str]:
        """Push entries, pickle self, submit jobs, collect failed entries."""
        queue = self._get_queue()
        self._clear_queues(queue)

        self._copy_entries_to_log_dir()
        self._enqueue_entries(queue)

        if queue.length(self.main_queue) == 0:
            raise BatchRunException(
                "No entries to process — queue is empty after enqueuing."
            )

        self._pickle_obj()
        self._log_cluster_config()

        actual_jobs = min(n_jobs, queue.length(self.main_queue))
        self._submit_to_farm(queue, actual_jobs, memory, self.main_queue)

        # Any entries still in main_queue were not picked up → treat as failed
        unprocessed = queue.get_all(self.main_queue)
        if unprocessed:
            logger.warning(
                f"{len(unprocessed)} entries were not picked up — "
                "moving to failed queue"
            )
            for item in unprocessed:
                queue.push(self.failed_queue, item)

        failed = queue.get_all(self.failed_queue)
        if failed:
            logger.info(f"{len(failed)} entries failed")
        return failed

    # ---------------------------------------------------------------------------
    # Cluster helpers
    # ---------------------------------------------------------------------------

    def process_remote(self) -> int:
        """Run on a cluster node: pop entries from the queue and process them.

        Called by :mod:`pdbe_sifts.base.process_remote` after unpickling self.

        Returns:
            Number of entries in the failed queue when this job exits.
        """
        from pdbe_sifts.base.queues.batchable_queue import BatchableQueue
        from pdbe_sifts.base.report import EntryStatusType, log_entry

        self.before_job_start()
        queue = BatchableQueue().get_queue()

        while True:
            entry = queue.pop_and_push(self.main_queue, self.failed_queue)
            if not entry:
                break

            try:
                self._process_entry_with_timeout(entry)
                queue.remove_item(self.failed_queue, entry)
                log_entry(entry, self._task_name, EntryStatusType.SUCCESS)
            except Exception:
                exc_type, exc_value, _ = sys.exc_info()
                log_entry(
                    entry,
                    self._task_name,
                    EntryStatusType.FAILED,
                    exc_type,
                    exc_value,
                )
                logger.error(f"Failed for entry {entry}:", exc_info=True)

        self.after_job_end()
        return queue.length(self.failed_queue)

    def _process_entry_with_timeout(self, entry_id: str) -> None:
        """Process an entry, respecting known_exceptions and entry_timeout_secs."""
        try:
            if not self.entry_timeout_secs or self.entry_timeout_secs <= 0:
                self.process_entry(entry_id)
            else:
                # Use a thread (not a subprocess) so the cluster node's DuckDB
                # connection (set up in before_job_start) is accessible.
                from multiprocessing.dummy import Process as ThreadProcess

                p = ThreadProcess(target=self.process_entry, args=(entry_id,))
                p.start()
                p.join(self.entry_timeout_secs)
                if p.is_alive():
                    p.join()  # thread cannot be killed — just wait briefly
                    raise TimeoutError(
                        f"Entry {entry_id!r} timed out after "
                        f"{self.entry_timeout_secs}s"
                    )
        except Exception:
            if entry_id not in self.known_exceptions:
                raise
            logger.warning(
                f"Entry {entry_id!r} is a known exception — failure ignored"
            )

    def _pickle_obj(self) -> None:
        pickle_dir = self._get_pickle_dir()
        pickle_dir.mkdir(parents=True, exist_ok=True)
        self.pickle_file = Path(pickle_dir, self.main_queue)
        with open(self.pickle_file, "wb") as fp:
            pickle.dump(self, fp)  # nosec
        logger.debug(f"Pickled self to {self.pickle_file}")

    def _submit_to_farm(self, queue, jobs: int, memory: int, job_name: str) -> None:
        """Submit *jobs* cluster jobs that each call process_remote.py."""
        from pdbe_sifts.base.process_remote import __file__ as process_remote_script

        remote_command = (
            f"{sys.executable} {process_remote_script} {self.pickle_file}"
        )
        logger.info(f"Submitting to cluster: {remote_command}")

        # Resolve cluster queue (may use GPU queue override)
        use_queue = self.cluster_queue
        if not use_queue:
            try:
                from pdbe_sifts.config import load_config
                use_queue = load_config().farm.cluster_queue or "normal"
            except Exception:
                use_queue = "normal"

        timeout_mins = None
        if self.entry_timeout_secs and self.entry_timeout_secs > 0:
            timeout_mins = self.entry_timeout_secs // 60

        self.executor.submit(
            job_name,
            remote_command,
            queue=use_queue,
            log_dir=self.log_dir,
            jobs=jobs,
            memory=memory,
            cores=self.cores,
            gpus=self.gpus,
            timeout_mins=timeout_mins,
        )

    def _clear_queues(self, queue) -> None:
        queue.delete(self.main_queue)
        queue.delete(self.failed_queue)

    def _enqueue_entries(self, queue) -> None:
        for entry in self.entries:
            queue.push(self.main_queue, entry)
        logger.info(
            f"{queue.length(self.main_queue)} entries pushed to queue "
            f"'{self.main_queue}'"
        )

    def _copy_entries_to_log_dir(self) -> None:
        ids_file = Path(self.log_dir, "all_ids.list")
        with open(ids_file, "w") as fh:
            for entry in self.entries:
                print(entry, file=fh)
        logger.info(f"Entry list written to {ids_file} ({len(self.entries)} entries)")

    def _log_cluster_config(self) -> None:
        logger.info("::params:: Cluster parameters")
        logger.info(f"  Jobs        : {self.jobs}")
        logger.info(f"  Cores       : {self.cores}")
        logger.info(f"  Init memory : {self.initial_memory} MB")
        logger.info(f"  Retry memory: {self.retry_memory} MB")
        logger.info(f"  Queue       : {self.cluster_queue}")
        logger.info(f"  Log dir     : {self.log_dir}")
        logger.info(f"  Timeout     : {self.entry_timeout_secs}s")
        logger.info(f"  Entries     : {len(self.entries)}")
        logger.info("::endparams::")

    def load_exceptions(self) -> None:
        """Load :attr:`exceptions_file` into :attr:`known_exceptions`."""
        if self.exceptions_file and Path(self.exceptions_file).exists():
            with open(self.exceptions_file) as fh:
                self.known_exceptions = {line.strip() for line in fh}
            logger.info(
                f"Loaded {len(self.known_exceptions)} known exceptions "
                f"from {self.exceptions_file}"
            )

    def write_failed_entries_to_file(self) -> None:
        """Write :attr:`failed_entries` to ``<log_dir>/failed.list``."""
        failed_file = Path(self.log_dir, "failed.list")
        with open(failed_file, "w") as fh:
            for entry in sorted(self.failed_entries):
                print(entry, file=fh)
        logger.info(f"Failed entries written to {failed_file}")

    def log_failed_entries(self) -> None:
        """Log the first few failed entries and grep logs for their errors."""
        count = min(5, len(self.failed_entries))
        logger.warning(
            f"{len(self.failed_entries)} entries failed. "
            f"First {count}: {list(self.failed_entries)[:count]}"
        )
        first = next(iter(self.failed_entries))
        grep_cmd = (
            f"grep -r {first} {self.log_dir}/*.out 2>/dev/null | grep -C 5 -i error"
        )
        try:
            result = subprocess.run(
                grep_cmd,
                shell=True,
                capture_output=True,
                check=True,
            )
            logger.info(
                f"Log excerpt for {first}:\n{result.stdout.decode().strip()}"
            )
        except subprocess.CalledProcessError:
            logger.debug(f"grep found no error lines in logs for {first}")

    def _send_summary(self) -> None:
        from pdbe_sifts.base.report import EntryStatusType, send_summary_message

        no_failed = len(self.failed_entries)
        unknown = self.failed_entries - self.known_exceptions
        status = (
            EntryStatusType.SUCCESS
            if no_failed == 0 or self.ignore_batch_failures
            else EntryStatusType.FAILED
        )
        send_summary_message(
            module=self._task_name,
            status=status,
            processed=len(self.entries),
            failed=len(unknown),
            known_exceptions=len(self.known_exceptions),
            log_dir=self.log_dir,
        )

    # ---------------------------------------------------------------------------
    # mmCIF revision utility (ported from ORC Batchable)
    # ---------------------------------------------------------------------------

    def no_used_cif_category_modified(self, cif_file: str) -> bool:
        """Return ``True`` if the CIF file has no modifications to monitored categories.

        Useful for skipping entries whose used mmCIF categories have not been
        updated since the last release.

        Scenarios:

        1. :attr:`used_cif_categories` is empty → return ``False`` (process it).
        2. Missing ``_pdbx_audit_revision_history`` block → return ``False``.
        3. No future-dated revision → return ``False``.
        4. Future revision exists but none of the modified categories are in
           :attr:`used_cif_categories` → return ``True`` (skip).

        Args:
            cif_file: Path to the mmCIF file.

        Returns:
            ``True`` if the entry can safely be skipped.
        """
        import gemmi

        if not self.used_cif_categories:
            logger.debug("no_used_cif_category_modified: no categories configured")
            return False

        used = {cat.lstrip("_") for cat in self.used_cif_categories}

        block = gemmi.cif.read(str(cif_file)).sole_block()
        history = block.find(
            "_pdbx_audit_revision_history.", ["ordinal", "revision_date"]
        )
        future_ordinals = [
            row["ordinal"]
            for row in history
            if self._is_future_date(row["revision_date"])
        ]

        if not future_ordinals:
            logger.info("no_used_cif_category_modified: no future revisions")
            return False

        categories = block.find(
            "_pdbx_audit_revision_category.", ["revision_ordinal", "category"]
        )
        if not categories:
            logger.info(
                "no_used_cif_category_modified: "
                "_pdbx_audit_revision_category not found"
            )
            return False

        modified = {
            row["category"]
            for row in categories
            if row["revision_ordinal"] in future_ordinals
        }
        overlap = modified & used
        if not overlap:
            logger.debug(f"Modified categories {modified} do not overlap with {used}")
            logger.info("no_used_cif_category_modified: none of the modified "
                        "categories are used — skipping")
            return True

        logger.info(f"no_used_cif_category_modified: overlap found: {overlap}")
        return False

    @staticmethod
    def _is_future_date(date_str: str) -> bool:
        try:
            return datetime.strptime(date_str, "%Y-%m-%d") > datetime.now()
        except ValueError:
            return False

    # ---------------------------------------------------------------------------
    # load_entries
    # ---------------------------------------------------------------------------

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

    # ---------------------------------------------------------------------------
    # CLI entry point
    # ---------------------------------------------------------------------------

    def add_arguments(self, parser: argparse.ArgumentParser) -> None:
        """Override to add custom CLI arguments.

        Called with the **root** ArgumentParser before subparsers are added,
        so arguments defined here are available in both ``single`` and
        ``batch`` subcommands.
        """

    def main(self, argv: Optional[list[str]] = None) -> None:
        """Parse CLI arguments and dispatch to the appropriate run method.

        Subcommands::

            single  --entry <id>
            batch   --list <file_or_csv> [--workers N] [--no-retry]
                    [--failure-threshold F] [--timeout S]
                    [--run-local]
                    [--jobs J] [--memory M] [--retry-memory R]
                    [--cores C] [--cluster lsf|slurm]
                    [--cluster-queue Q] [--queue-name N] [--log-dir D]

        Args:
            argv: Argument list to parse.  Uses ``sys.argv[1:]`` when ``None``.
        """
        root = argparse.ArgumentParser(description=self.__class__.__name__)
        self.add_arguments(root)
        sub = root.add_subparsers(dest="mode", required=True)

        # --- single subcommand ---
        single_p = sub.add_parser("single", help="Process one entry.")
        single_p.add_argument(
            "--entry", "-e", required=True, help="Entry ID to process."
        )

        # --- batch subcommand ---
        batch_p = sub.add_parser(
            "batch", help="Process multiple entries (local or cluster)."
        )
        batch_p.add_argument(
            "--list", "-l",
            dest="entry_list",
            help=(
                "Path to a file of entry IDs (one per line) or a comma-separated "
                "list.  Uses entry_file_path if not provided."
            ),
        )

        # Local parallel options
        local_grp = batch_p.add_argument_group("Local execution")
        local_grp.add_argument(
            "--run-local",
            action="store_true",
            default=True,
            help="Use local multiprocessing Pool (default).  "
                 "Ignored when --cluster is set.",
        )
        local_grp.add_argument(
            "--workers", "-w",
            type=int,
            default=self.workers,
            help=f"Worker processes for local mode (default: {self.workers}).",
        )
        local_grp.add_argument(
            "--no-retry",
            action="store_true",
            help="Disable automatic retry for failed entries.",
        )
        local_grp.add_argument(
            "--failure-threshold",
            type=float,
            default=self.failure_threshold,
            dest="failure_threshold",
            help=(
                "Max allowed failure ratio (<1) or count (>=1). "
                "Default: %(default)s."
            ),
        )
        local_grp.add_argument(
            "--timeout",
            type=int,
            default=self.entry_timeout_secs,
            dest="entry_timeout_secs",
            help="Per-entry timeout in seconds (local mode).  No timeout if not set.",
        )

        # Cluster options
        cluster_grp = batch_p.add_argument_group("Cluster execution")
        cluster_grp.add_argument(
            "--cluster",
            choices=["lsf", "slurm"],
            default=None,
            metavar="TYPE",
            help="Submit to cluster (lsf or slurm).  Activates cluster mode.",
        )
        cluster_grp.add_argument(
            "--jobs", "-j",
            type=int,
            default=self.jobs,
            help=f"Number of parallel cluster jobs (default: {self.jobs}).",
        )
        cluster_grp.add_argument(
            "--memory",
            type=int,
            default=self.initial_memory,
            dest="memory",
            help=f"Initial memory per cluster job in MB (default: {self.initial_memory}).",
        )
        cluster_grp.add_argument(
            "--retry-memory",
            type=int,
            default=self.retry_memory,
            dest="retry_memory",
            help=f"Memory for retry jobs in MB (default: {self.retry_memory}).",
        )
        cluster_grp.add_argument(
            "--cores",
            type=int,
            default=self.cores,
            help=f"CPU cores per cluster job (default: {self.cores}).",
        )
        cluster_grp.add_argument(
            "--cluster-queue",
            default=self.cluster_queue or None,
            dest="cluster_queue",
            help="Cluster queue/partition name.",
        )
        cluster_grp.add_argument(
            "--queue-name",
            default=None,
            dest="queue_name",
            help="Messaging queue name (auto-generated if not given).",
        )
        cluster_grp.add_argument(
            "--log-dir",
            default=None,
            dest="log_dir",
            help="Directory for cluster job stdout/stderr.",
        )

        args = root.parse_args(argv)
        logger.info(vars(args))

        if args.mode == "single":
            self.run_single(args.entry)
            return

        # batch mode
        self.failure_threshold = args.failure_threshold
        self.entry_timeout_secs = args.entry_timeout_secs

        entries = (
            _parse_entry_list(args.entry_list)
            if args.entry_list
            else self.load_entries()
        )

        use_cluster = args.cluster is not None
        if use_cluster:
            self.run_cluster_batch(
                entries,
                jobs=args.jobs,
                memory=args.memory,
                retry_memory=args.retry_memory,
                cluster=args.cluster,
                cluster_queue=args.cluster_queue,
                queue_name=args.queue_name,
                log_dir=args.log_dir,
                no_retry=args.no_retry,
            )
        else:
            self.run_batch(
                entries,
                workers=args.workers,
                retry=not args.no_retry,
            )

    # ---------------------------------------------------------------------------
    # Internal helpers (local mode)
    # ---------------------------------------------------------------------------

    def _run_parallel(self, entries: list[str], n_workers: int) -> list[str]:
        """Distribute entries across worker processes; return list of failed IDs.

        Uses multiprocessing.Pool with imap_unordered so that:
        - Each entry is an independent task (one crash → one entry lost, not a chunk)
        - Crashed workers are automatically restarted by the Pool
        - Workers are recycled after SIFTS_MAXTASKS_PER_CHILD entries to release memory
        """
        failed: list[str] = []
        completed: set[str] = set()
        n_workers = max(1, n_workers)
        maxtasks = int(os.environ.get("SIFTS_MAXTASKS_PER_CHILD", 500)) or None

        with Pool(
            processes=n_workers,
            initializer=_init_worker,
            initargs=(self,),
            maxtasksperchild=maxtasks,
        ) as pool:
            results = pool.imap_unordered(_run_entry, entries, chunksize=1)
            while True:
                try:
                    entry_id, exc_tb = next(results)
                    completed.add(entry_id)
                    if exc_tb:
                        logger.error(f"Failed [{entry_id}]:\n{exc_tb}")
                        failed.append(entry_id)
                except StopIteration:
                    break
                except Exception as exc:
                    # A worker was killed (OOM/SIGKILL). The Pool has already
                    # restarted a replacement worker. We cannot know which
                    # entry_id was being processed, so we log and continue
                    # consuming results from the still-live iterator.
                    logger.error(
                        f"Worker crashed ({type(exc).__name__}: {exc}). "
                        "Pool has restarted a replacement worker.",
                        exc_info=True,
                    )

        # Entries that never returned a result were on the crashed worker
        crashed = set(entries) - completed
        if crashed:
            logger.warning(
                f"{len(crashed)} entries lost to worker crash(es): "
                f"{sorted(crashed)[:5]}{'...' if len(crashed) > 5 else ''}"
            )
            failed.extend(sorted(crashed))

        return failed

    def _check_failures(self, failed: list[str], total: int) -> None:
        if not failed:
            logger.info("All entries processed successfully.")
            return
        if self.ignore_batch_failures or self.failure_threshold <= 0:
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
        if self.failed_entries_file is not None:
            path = Path(self.failed_entries_file)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("\n".join(failed) + "\n")
            logger.info(f"Failed entries written to {path}")
