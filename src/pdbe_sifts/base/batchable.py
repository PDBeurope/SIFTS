import argparse
import os
import pickle  # nosec
import subprocess
import sys
from abc import ABC, abstractmethod
from datetime import datetime
from enum import Enum
from multiprocessing.dummy import Process
from pathlib import Path
from typing import Iterable

import gemmi
from pdbe_sifts.config import load_config

from pdbe_sifts.base.exceptions import (
    BatchRunException,
    EntryTimedOutException,
    TooManyFailedEntries,
)
from pdbe_sifts.base.executors import get_executor
from pdbe_sifts.base.log import logger
from pdbe_sifts.base.parser import CORES, INIT_MEM, JOBS, RETRY_MEM, GPUS
from pdbe_sifts.base.process_remote import __file__ as process_remote_script
from pdbe_sifts.base.queues.batchable_queue import BatchableQueue
from pdbe_sifts.base.queues.iqueue import IQueue
# from pdbe_sifts.base.report import EntryStatusType, log_entry, send_summary_message


class Modes(Enum):
    SINGLE = "SINGLE"
    BATCH = "BATCH"


conf = load_config()


class Batchable(ABC):
    """Abstract Base Class for tasks that can be parallelized.

    the parallelization can be either locally or on remote cluster (LSF for now).
    The subclass can additionally set the attributes below (shown with defaults):

    Raises:
        ValueError: If neither --list option specified nor self.entry_file_path
            set to a valid file, in batch mode
        BatchRunException: If some entries did not process successfully.
            Overridden by ignore_batch_failures and failure_threshold.
        TooManyFailedEntries: If more than `failure_threshold` entries failed
    """  # noqa: E501,B950

    initial_memory: int = INIT_MEM
    """Initial memory in MB to submit with."""

    retry_memory: int = RETRY_MEM
    """Memory in MB to retry failed entries with."""

    jobs: int = JOBS
    """Jobs/threads to submit in parallel"""

    cores: int = CORES
    """Number of cores to request to use/request for each job"""

    gpus: str = GPUS
    """Number of gpus to request to use/request for each job"""

    cluster_queue: str = conf.farm.cluster_queue
    """Name of cluster queue to submit to."""

    entries: list[str] = []
    """List of entries to process."""

    failed_entries: set = set()
    """Set of failed entries"""

    entry_file_path: str | None = None
    """File with list of entries to process. Optional"""

    no_retry: bool = False
    """If True, do no retry if there are failures in batch run"""

    ignore_batch_failures: bool = False
    """If True, does not raise exception if batch run fails. Overrides failure_threshold."""

    known_exceptions: set = set()
    """List of known exceptions. Optional"""

    exceptions_file: str | None = None
    """File with list of entries of known exceptions. Optional"""

    log_dir: str = ""
    """Directory to write logs to. Defaults to conf.location.work.bsub_logs_dir."""

    mode: Modes = Modes.SINGLE
    """Mode to run in. SINGLE or BATCH"""

    filter: range | None = None  # noqa: A003
    """Range of entries to process. Optional"""

    failure_threshold: float = 0.0
    """Number or fraction of failures to tolerate before raising exception."""

    entry_timeout_secs = None
    """Number of seconds to wait for entry to be processed before timing out."""

    pickle_dir = Path(
        os.getenv("SIFTS_PICKLE_DIR") or conf.location.work.pickle_dir or Path.home()
    )
    """Directory to store pickled objects. Defaults to envvar SIFTS_PICKLE_DIR.

    If not set, defaults to user's home directory.
    """

    pickle_file = None
    """File to pickle self to. Optional. defaults to <pickle_dir>/<main_queue>"""

    used_cif_categories: Iterable[str] = set()

    executor = None

    _start_timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

    _task_name = "task"

    def _log_config(self):
        logger.info("::params::Parameters:")
        logger.info(f"Cores: {self.cores}")
        logger.info(f"Gpus: {self.gpus}")
        logger.info(f"Initial Memory: {self.initial_memory}")
        logger.info(f"Retry Memory: {self.retry_memory}")
        logger.info(f"Jobs: {self.jobs}")
        logger.info(f"Cluster Queue: {self.cluster_queue}")
        logger.info(f"Log Directory: {self.log_dir}")
        logger.info(f"Filter: {self.filter}")
        logger.info(f"Failure Threshold: {self.failure_threshold}")
        logger.info(f"Entry Timeout: {self.entry_timeout_secs}")
        logger.info(f"Run Local: {self.run_local}")
        logger.info(f"List of Entries: {self.entry_file_path}")
        logger.info(f"Exceptions File: {self.exceptions_file}")
        logger.info(f"No of Known Exceptions: len({self.known_exceptions})")
        logger.info("::endparams::")

    def _clear_queues(self, queue: IQueue):
        """Clears failed and main queues if exist."""
        queue.delete(self.main_queue)
        queue.delete(self.failed_queue)

    def load_entries(self):
        """Loads self.entries list from a file as specified by `entries_file_path`.

        Override this to have a custom way of populating self.entries.
        This does not load into the queue. Use `_enqueue_entries` for that.
        """
        if self.entry_file_path:
            if not Path(self.entry_file_path).exists():
                raise FileNotFoundError(
                    f"File {self.entry_file_path} must exists and be readable"
                )
            if not Path(self.entry_file_path).lstat().st_size > 0:
                raise ValueError(
                    f"List of entries file {self.entry_file_path} cannot be empty"
                )

            with open(self.entry_file_path) as f:
                self.entries = [line.strip() for line in f]

        if not self.entries:
            raise ValueError("No entries to process.")
        logger.info(f"{len(self.entries)} items added for processing")

    def _copy_entries_to_log_dir(self):
        ids_copy_file = os.path.join(self.log_dir, "all_ids.list")
        with open(ids_copy_file, "w+") as f:
            logger.info(f"Will process {len(self.entries)} entries")
            for entry in self.entries:
                print(entry, file=f)
        logger.info(f"Ids to process copied to: {ids_copy_file}")

    def _enqueue_entries(self, queue: IQueue):
        """Enqueues entries to messaging queue."""

        for e in self.entries:
            queue.push(self.main_queue, e)

        count_items = queue.length(self.main_queue)
        logger.info(f"{count_items} items added to queue {self.main_queue}")

    def submit_batch(self):
        """Submit batch jobs the compute cluster.

        Name of job in cluster is same as self.main_queue
        """
        self._copy_entries_to_log_dir()

        self.failed_queue = f"{self.main_queue}-failed"

        queue = BatchableQueue().get_queue()

        self._clear_queues(queue)

        logger.debug(f"Loading {len(self.entries)} entries into messaging queue")
        self._enqueue_entries(queue)

        if not self.entries:
            raise BatchRunException(
                "No entries to process. What are the chances? so failing"
            )

        self._pickle_obj()

        self._log_config()

        try:
            jobs = min(self.jobs, queue.length(self.main_queue))

            self.failed_entries = set(
                self._submit_to_farm(queue, jobs, self.initial_memory, self.main_queue)
            )

            if not self.no_retry and self.failed_entries and self.before_retry():
                logger.info(f"Retrying {len(self.failed_entries)} failed entries")

                queue.rename(self.failed_queue, self.main_queue)

                jobs = min(len(self.failed_entries), self.jobs)

                self.failed_entries = self._submit_to_farm(
                    queue, jobs, self.retry_memory, f"{self.main_queue}-failed"
                )

            if self.failed_entries:
                self.write_failed_entries_to_file()
                self.log_failed_entries()

            else:
                logger.info("Processed all entries successfully")

        finally:
            logger.debug("removing  pickle file")
            if self.pickle_file and os.path.exists(self.pickle_file):
                os.remove(self.pickle_file)

    def log_failed_entries(self):
        count = min(5, len(self.failed_entries))
        logger.info(
            f"First {count} failed entries: {list(self.failed_entries)[:count]}"
        )
        one_failed = list(self.failed_entries)[0]
        # Log the first failed entry in detail
        logger.info(f"Logging details for failed entry {one_failed}")
        # grep for the entry in all log files in the log_dir and extract the error lines
        # with context
        grep_cmd = f"grep -r {one_failed} {self.log_dir}/*.out | grep -C 10 -i error"
        logger.info(f"Running: {grep_cmd}")
        try:
            grep_result = subprocess.run(
                grep_cmd,
                shell=True,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            logger.info(
                f"Failed log for {one_failed} : {grep_result.stdout.decode().strip()}"
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"Grep command failed with error: {e.stderr.decode().strip()}")
            print(e.stdout.decode().strip())

    def write_failed_entries_to_file(self):
        failed_file = Path(self.log_dir, "failed.list")
        with open(failed_file, "w+") as f:
            for e in self.failed_entries:
                print(e, file=f)
        logger.info(f"List of failed entries written to {failed_file}")

    def job_cleanup(self):  # noqa: B027
        """Override this to do any job clean up"""
        pass

    def _pickle_obj(self):
        self.pickle_file = self.pickle_file or Path(self.pickle_dir, self.main_queue)
        with open(self.pickle_file, "wb") as fp:
            logger.debug("Dumping pickle")
            pickle.dump(self, fp)

    def _submit_to_farm(self, queue: IQueue, jobs: int, memory: int, job_name: str):
        """Sends a job array to the farm which uses the pickled object to run from queue."""
        remote_command = f"{sys.executable} {process_remote_script} {self.pickle_file}"
        logger.info(f"Submitting command to farm: {remote_command}")

        cluster_queue = self.cluster_queue or conf.farm.cluster_queue
        timeout_mins = (
            self.entry_timeout_secs // 60
            if self.entry_timeout_secs and self.entry_timeout_secs > 0
            else None
        )

        gpus = "0"
        if self.gpus and not str(self.gpus) == "0":
            gpus = self.gpus
            cluster_queue = "short"

        self.executor.submit(
            job_name,
            remote_command,
            queue=cluster_queue,
            log_dir=str(self.log_dir),
            jobs=jobs,
            memory=memory,
            cores=self.cores,
            gpus=gpus,
            timeout_mins=timeout_mins,
        )

        unprocessed = queue.get_all(self.main_queue)
        if unprocessed:
            logger.warning(
                f"Some {len(unprocessed)} entries are unprocessed. "
                "Merging them to failed entries"
            )
            for item in unprocessed:
                queue.push(self.failed_queue, item)

        failed_entries = queue.get_all(self.failed_queue)
        if failed_entries:
            logger.info(f"{len(failed_entries)} failed entries")

        return failed_entries

    def pre_run(self):  # noqa: B027
        """Override this to do setup before run is called."""
        pass

    def before_job_start(self):  # noqa: B027
        """Override this to specify per job pre-process logic.

        This is executed before the start of each job on the farm in batch mode,
        or before process_entry for single mode. Use this to do setup of values
        that should be used across the job e.g. a database connection that will
        be used for all entries processed by this job or, a working area for that
        is specific to the job.
        """
        pass

    @abstractmethod
    def process_entry(self, entry_id: str):
        """Abstract method to process single entry.

        Each subclass should implement this method to process a single entry.

        Args:
            entry_id: ID to process e.g. CCD ID, PDB ID, PRD ID, Uniprot Accession
        """
        pass

    def after_job_end(self):  # noqa: B027
        """Override this method to specify per job post-process logic.

        This is executed after all entries have been executed for a specific farm job.
        for single mode, it is executed after `process_entry`. Use this for cleanup of
        temp files, database connections etc specified by before_job_start.
        """
        pass

    def post_run(self):  # noqa: B027
        """Override this to do summary and cleanup after run completes."""
        pass

    def before_retry(self) -> bool:
        """Override this method to add logic to be run before retry.

        Use this to conditionally interrupt retries to either set some specific values
        or cancel the retry.

        Returns:
            Return False to cancel retrying, True to go ahead with retry.
        """
        return True

    def process_single(self, entry_id: str):
        """Process entry in single mode.

        Wrapper to run pre- and post-process tasks.

        Args:
            entry_id: ID to process passed to `process_entry`
        """
        self.before_job_start()
        self.__process_entry(entry_id)
        self.after_job_end()

    def process_remote(self) -> int:
        """Process entries, picking them from a queue."""
        self.before_job_start()

        queue: IQueue = BatchableQueue().get_queue()
        while True:
            entry = queue.pop_and_push(self.main_queue, self.failed_queue)
            if not entry:
                break

            try:
                self.__process_entry(entry)
                queue.remove_item(self.failed_queue, entry)
                # log_entry(entry, self._task_name, EntryStatusType.SUCCESS)
            except Exception:
                exc_type, exc_value, _ = sys.exc_info()
                # log_entry(
                #     entry,
                #     self._task_name,
                #     EntryStatusType.FAILED,
                #     exc_type,
                #     exc_value,
                # )

                logger.error(f"Failed for entry {entry}:", exc_info=True)

        self.after_job_end()

        return queue.length(self.failed_queue)

    def raise_(self, ex):
        raise ex

    def __process_entry(self, entry_id):
        try:
            if not self.entry_timeout_secs or self.entry_timeout_secs <= 0:
                self.process_entry(entry_id)
            else:
                # Actually, its a thread, not a process. (ref: import multiprocessing.dummy)
                p = Process(
                    target=self.process_entry,
                    args=(entry_id,),
                )

                # Inject sub-method to stop the thread
                p.stop = lambda: self.raise_(
                    EntryTimedOutException(f"Timeout for {entry_id}")
                )
                p.start()
                p.join(self.entry_timeout_secs)

                if p.is_alive():
                    p.stop()
                    p.join()

        except Exception:
            if entry_id not in self.known_exceptions:
                raise
            logger.warning(
                f"Entry {entry_id} is in known_exceptions. Will ignore failure"
            )

    def local_batch(self):
        """Process the batch locally when --run-local is specified.

        Uses multiprocessing if cores > 1. Falls back to mutithreading otherwise

        Returns:
            list -- results of each entry processing as a list
        """
        self._copy_entries_to_log_dir()

        if self.cores > 1:
            from multiprocessing import Pool

            self.jobs = self.cores
        else:
            from multiprocessing.dummy import Pool

        p = Pool(self.jobs)
        results = p.map(self.process_single, self.entries)
        return results

    def load_exceptions(self):
        if self.exceptions_file and Path(self.exceptions_file).exists():
            with open(self.exceptions_file) as f:
                self.known_exceptions = {line.strip() for line in f}

    def setup(self, args):
        """Hydrates the object members from passed command line arguments."""

        # Convert argparse.Namespace to Dictionary
        if isinstance(args, argparse.Namespace):
            args = vars(args)

        self._task_name = self.__class__.__name__

        self.mode = Modes(args["subcommand"].upper())
        if self.mode == Modes.SINGLE:
            self.entries = [args["entry"]]
            return
        self.filter = args["filter"]

        self.ignore_batch_failures = (
            args["ignore_batch_failures"] or self.ignore_batch_failures
        )
        self.failure_threshold = float(
            args["failure_threshold"] or self.failure_threshold
        )
        self.entry_timeout_secs = (
            int(args["entry_timeout_secs"]) if args["entry_timeout_secs"] else None
        )

        self.executor = get_executor(args.get("executor"))
        self.initial_memory = args["initial_memory"] or self.initial_memory
        self.cores = args["cores"] or self.cores
        self.gpus = args["gpus"] or self.gpus or "0"
        self.retry_memory = args["retry_memory"] or self.retry_memory
        self.retry_memory = max(self.initial_memory, self.retry_memory)
        self.no_retry = args["no_retry"]
        self.jobs = args["no_of_jobs"] or self.jobs
        self.cluster_queue = args["cluster_queue"] or self.cluster_queue
        self.main_queue = args["main_queue_name"] or "_".join(
            [self._task_name.lower(), self._start_timestamp]
        )
        log_base_dir = (
            args["log_dir"]  # Specified by user
            or os.getenv("SIFTS_EXECUTOR_LOG_DIR")  # Specified in env
            or self.log_dir  # Specified in subclass
            or conf.executor.log_base_dir  # Default
        )
        self.log_dir = (
            f"{Path(log_base_dir, self._task_name.lower(), self._start_timestamp)}"
        )

        self.pickle_dir.mkdir(parents=True, exist_ok=True)

        if args["list"] and isinstance(args["list"], str):
            self.entry_file_path = args["list"]
            self.entries = []
        if args["list"] and isinstance(args["list"], list):
            self.entry_file_path = None
            self.entries = args["list"]
        self.run_local = True if args["run_local"] else False

    def _is_future_date(self, date: str) -> bool:
        """Check if the date is in the future.

        Args:
            date (str): Date in the format "YYYY-MM-DD"

        Returns:
            bool: True if the date is in the future, False otherwise.
        """
        return datetime.strptime(date, "%Y-%m-%d") > datetime.now()

    def no_used_cif_category_modified(self, cif_file: str) -> bool:
        """Returns true if the CIF file has no modifications in the categories used by the task.

        Scenarios
        1. Missing _pdbx_audit_revision_history or _pdbx_audit_revision_category, return False
        2. `self.used_cif_categories` is empty, return False.
        3. Latest revision is not in the future, return False.
        4. Modified categories exist, but none of them are used by task, return True.

        Args:
            cif_file (str): Path to the CIF file

        Returns:
            bool: True if the CIF file has no modifications in the categories used by the task.
                  False otherwise including if the categories are not found in the CIF file.
        """
        if not self.used_cif_categories:
            logger.debug("No categories to check for modifications")
            return False

        # Remove leading _ in self.used_cif_categories if present
        # since the categories in pdbx_ausit_category are without leading _
        self.used_cif_categories = {cat.lstrip("_") for cat in self.used_cif_categories}

        block = gemmi.cif.read(str(cif_file)).sole_block()
        history = block.find(
            "_pdbx_audit_revision_history.", ["ordinal", "revision_date"]
        )
        ordinals = [
            row["ordinal"]
            for row in history
            if self._is_future_date(row["revision_date"])
        ]

        if not ordinals:
            logger.info("No future revisions found")
            return False

        categories = block.find(
            "_pdbx_audit_revision_category.", ["revision_ordinal", "category"]
        )
        if not categories:
            logger.info("No pdbx_audit_revision_category found")
            return False

        modified_categories = {
            row["category"] for row in categories if row["revision_ordinal"] in ordinals
        }
        modified_used_categories = modified_categories.intersection(
            self.used_cif_categories
        )
        if not modified_used_categories:
            logger.debug(f"Modified categories: {', '.join(modified_categories)}")
            logger.debug(f"Used categories: {', '.join(self.used_cif_categories)}")
            logger.info("None of the modified categories are used.")
            return True

        logger.info(f"Modified categories found: {', '.join(modified_used_categories)}")
        return False

    def run(self, args: argparse.Namespace | dict):
        """Run application based on args options.

        If dictionary is specified, the keys should be the same as the base parser destinations.

        Args:
            args: Parsed Args in either dictionary

        Raises:
            KeyError: if dictionary inoput provided is missing required key
        """
        self.setup(args)

        self.pre_run()
        self.load_exceptions()

        try:
            if self.mode == Modes.SINGLE:
                for entry in self.entries:
                    self.process_single(entry)

            elif self.mode == Modes.BATCH:
                os.makedirs(self.log_dir, exist_ok=True)
                self.load_entries()
                if self.run_local:
                    self.local_batch()
                else:
                    self.submit_batch()
                self.check_failure()

        finally:
            self.post_run()
            if self.mode == Modes.BATCH:
                pass
                # self.send_summary()

    def check_failure(self):
        if not self.failed_entries:
            logger.debug("All entries processed successfully")
            return

        failures = len(self.failed_entries)
        if self.ignore_batch_failures:
            logger.warning(f"Batch failed with {failures} entries. Ignoring failures")
            return

        allowed_failures = (
            int(self.failure_threshold * len(self.entries))
            if self.failure_threshold < 1
            else self.failure_threshold
        )

        logger.debug(f"Threshold: {self.failure_threshold} ==> {allowed_failures}")
        logger.debug(f"Failed: {failures}")

        if failures > allowed_failures:
            msg = f"failures ({failures}) > threshold({allowed_failures})"
            raise TooManyFailedEntries(msg)

    # def send_summary(self):
    #     no_of_entries = len(self.entries)
    #     no_of_known_exceptions = len(self.known_exceptions)
    #     no_failed = len(self.failed_entries)
    #     unknown_exceptions = set(self.failed_entries).difference(self.known_exceptions)

    #     status = EntryStatusType.SUCCESS
    #     if no_failed > 0 and not self.ignore_batch_failures:
    #         status = EntryStatusType.FAILED

    #     send_summary_message(
    #         self._task_name,
    #         status,
    #         no_of_entries,
    #         len(unknown_exceptions),
    #         no_of_known_exceptions,
    #         self.log_dir,
    #     )
