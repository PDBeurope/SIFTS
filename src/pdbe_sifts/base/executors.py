"""Cluster job executors for LSF and Slurm.

All config values are read **lazily** (inside methods, not at import time) so
that this module can be imported on machines without a config file.

The public API is :func:`get_executor`, which returns an :class:`Executor`
subclass based on the ``executor.type`` config key or the ``--cluster`` CLI
flag.
"""

import os
import subprocess
import time
from abc import ABC, abstractmethod
from collections import defaultdict
from enum import Enum
from textwrap import dedent

from tqdm import tqdm

from pdbe_sifts.base.log import logger


class JobStatus(Enum):
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"
    OTHER = "OTHER"


class Executor(ABC):
    """Abstract executor interface."""

    job_ids: set[int]

    @abstractmethod
    def submit(self, job_name: str, command: str, **kwargs) -> None:
        """Submit *jobs* copies of *command* to the cluster and monitor until done."""

    @abstractmethod
    def get_job_status_by_id(self, job_id: int) -> str:
        """Return the raw status string for *job_id*."""

    @abstractmethod
    def map_job_status(self, status_text: str) -> JobStatus:
        """Map a raw status string to a :class:`JobStatus` enum value."""

    @abstractmethod
    def requeue_job(self, job_id: int) -> None:
        """Requeue a failed job."""

    def monitor(self, frequency: int = 5) -> None:
        """Block until all submitted jobs have finished, polling every *frequency* seconds."""
        with tqdm(total=len(self.job_ids), desc="Job completion", disable=None) as pbar:
            total_jobs = len(self.job_ids)
            completed_jobs = 0
            cancelled_jobs = 0
            failed_jobs: set[int] = set()

            while True:
                job_statuses: dict[JobStatus, set] = defaultdict(set)

                for job_id in list(self.job_ids):
                    status_text = self.get_job_status_by_id(job_id)
                    status = self.map_job_status(status_text)

                    if status == JobStatus.CANCELLED:
                        cancelled_jobs += 1
                        logger.warning(f"Job {job_id} was cancelled — not requeuing")
                        self.job_ids.discard(job_id)

                    elif status == JobStatus.COMPLETED:
                        completed_jobs += 1
                        self.job_ids.discard(job_id)
                        failed_jobs.discard(job_id)

                    elif status == JobStatus.FAILED:
                        try:
                            logger.debug(f"Job {job_id} failed — requeuing")
                            self.requeue_job(job_id)
                        except Exception:
                            logger.error(
                                f"Error requeuing job {job_id}", exc_info=True
                            )
                            failed_jobs.add(job_id)
                    else:
                        logger.debug(f"Job {job_id} status: {status}")
                        job_statuses[status].add(job_id)

                total_done = completed_jobs + cancelled_jobs + len(failed_jobs)
                pbar.update(total_done - pbar.n)

                if total_done == total_jobs:
                    logger.info("All jobs completed")
                    break

                time.sleep(frequency)


# ---------------------------------------------------------------------------
# LSF executor
# ---------------------------------------------------------------------------

class LSFJobManager(Executor):
    """Submit and monitor jobs on an LSF (bsub/bjobs) cluster."""

    def __init__(self):
        self.job_ids: set[int] = set()

    def submit(self, job_name: str, command: str, **kwargs) -> None:
        """Submit *jobs* independent copies of *command* via ``bsub`` and monitor."""
        cmd = self.build_command(job_name, command, **kwargs)
        logger.info(f"Job name: {job_name}")
        logger.info(f"Job options: {kwargs}")
        logger.debug(f"bsub command: {cmd}")

        for _ in range(kwargs.get("jobs", 1)):
            result = subprocess.run(cmd, check=True, capture_output=True)
            # bsub output: "Job <12345> is submitted to queue <normal>."
            job_id = int(result.stdout.decode().split()[1][1:-1])
            logger.debug(f"Submitted LSF job: {job_id}")
            self.job_ids.add(job_id)

        self.monitor(kwargs.get("monitor_frequency", 5))

    def build_command(self, job_name: str, command: str, **kwargs) -> list[str]:
        """Build a ``bsub`` command list from *kwargs*.

        Config defaults are read lazily so the method is safe to call even
        without a running config at import time.
        """
        from pdbe_sifts.config import load_config
        cfg = load_config()

        log_dir = (
            kwargs.get("log_dir")
            or os.getenv("SIFTS_EXECUTOR_LOG_DIR")
            or cfg.executor.log_base_dir
        )
        os.makedirs(log_dir, exist_ok=True)

        return [
            "bsub",
            "-J", job_name,
            "-q", str(kwargs.get("queue") or cfg.executor.default_queue),
            "-n", str(kwargs.get("cores") or cfg.executor.default_cores),
            "-M", str(kwargs.get("memory") or cfg.executor.default_memory),
            "-W", str(kwargs.get("timeout_mins") or cfg.executor.default_timeout),
            "-cwd", log_dir,
            "-o", kwargs.get("log_file", "lsf-%J.out"),
            command,
        ]

    def requeue_job(self, job_id: int) -> None:
        subprocess.run(["brequeue", "-e", str(job_id)], check=True)
        logger.debug(f"Re-queued LSF job {job_id}")

    def get_job_status_by_id(self, job_id: int) -> str:
        result = subprocess.run(
            ["bjobs", "-o", "stat", "-noheader", str(job_id)],
            check=True,
            capture_output=True,
        )
        status_text = result.stdout.decode().strip()
        if status_text == "EXIT":
            # Distinguish user-cancelled from real failures
            detail = subprocess.run(
                ["bjobs", "-l", str(job_id)], check=True, capture_output=True
            )
            if "TERM_OWNER" in detail.stdout.decode():
                status_text = "CANCELLED"
        return status_text

    def map_job_status(self, status_text: str) -> JobStatus:
        mapping = {
            "DONE": JobStatus.COMPLETED,
            "RUN": JobStatus.RUNNING,
            "PEND": JobStatus.RUNNING,
            "PSUSP": JobStatus.RUNNING,
            "EXIT": JobStatus.FAILED,
            "CANCELLED": JobStatus.CANCELLED,
        }
        return mapping.get(status_text, JobStatus.OTHER)


# ---------------------------------------------------------------------------
# Slurm executor
# ---------------------------------------------------------------------------

class SlurmJobManager(Executor):
    """Submit and monitor jobs on a Slurm (sbatch/squeue) cluster."""

    def __init__(self):
        self.job_ids: set[int] = set()

    def submit(self, job_name: str, command: str, **kwargs) -> None:
        """Submit *jobs* independent copies of *command* via ``sbatch`` and monitor."""
        cmd = self.build_command(job_name, command, **kwargs)
        logger.debug(f"sbatch command: {cmd}")

        for _ in range(kwargs.get("jobs", 1)):
            result = subprocess.run(cmd, check=True, capture_output=True)
            # sbatch output: "Submitted batch job 12345"
            job_id = int(result.stdout.decode().split()[-1])
            logger.debug(f"Submitted Slurm job: {job_id}")
            self.job_ids.add(job_id)

        self.monitor(kwargs.get("monitor_frequency", 5))

    def build_command(self, job_name: str, command: str, **kwargs) -> list[str]:
        """Build an ``sbatch`` command list from *kwargs*."""
        from pdbe_sifts.config import load_config
        cfg = load_config()

        log_dir = (
            kwargs.get("log_dir")
            or os.getenv("SIFTS_EXECUTOR_LOG_DIR")
            or cfg.executor.log_base_dir
        )
        os.makedirs(log_dir, exist_ok=True)

        # Write a tiny wrapper script for sbatch
        script_file = os.path.join(log_dir, f"{job_name}.sh")
        with open(script_file, "w") as f:
            f.write(dedent(f"""\
                #!/bin/bash
                set -e
                {command}
            """))

        return [
            "sbatch",
            "--job-name", job_name,
            "--partition", str(kwargs.get("queue") or cfg.executor.default_queue),
            "--cpus-per-task", str(kwargs.get("cores") or cfg.executor.default_cores),
            "--mem", str(kwargs.get("memory") or cfg.executor.default_memory),
            "--time", str(kwargs.get("timeout_mins") or cfg.executor.default_timeout),
            "--chdir", log_dir,
            "--open-mode", "append",
            "--output", kwargs.get("log_file", "slurm-%j.out"),
            script_file,
        ]

    def requeue_job(self, job_id: int) -> None:
        subprocess.run(["scontrol", "requeue", str(job_id)], check=True)
        logger.debug(f"Re-queued Slurm job {job_id}")

    def get_job_status_by_id(self, job_id: int) -> str:
        result = subprocess.run(
            ["squeue", "--noheader", "-t", "all", "--Format", "State",
             "--jobs", str(job_id)],
            check=True,
            capture_output=True,
        )
        return result.stdout.decode().strip()

    def map_job_status(self, status_text: str) -> JobStatus:
        if status_text in {"FAILED", "TIMEOUT", "OUT_OF_MEMORY"}:
            return JobStatus.FAILED
        elif status_text == "COMPLETED":
            return JobStatus.COMPLETED
        elif status_text in {"RUNNING", "PENDING", "CONFIGURING", "COMPLETING"}:
            return JobStatus.RUNNING
        elif status_text == "CANCELLED":
            return JobStatus.CANCELLED
        return JobStatus.OTHER


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

def get_executor(executor_type: str | None = None) -> Executor:
    """Return an :class:`Executor` for the given cluster type.

    The type is resolved in priority order:

    1. *executor_type* argument
    2. ``SIFTS_EXECUTOR`` environment variable
    3. ``executor.type`` in config.yaml

    Args:
        executor_type: ``"lsf"`` or ``"slurm"`` (case-insensitive).

    Raises:
        ValueError: If no executor type could be determined or the type is unknown.
    """
    if not executor_type:
        from pdbe_sifts.config import load_config
        cfg = load_config()
        executor_type = os.getenv("SIFTS_EXECUTOR") or cfg.executor.type

    if not executor_type:
        raise ValueError(
            "No cluster executor configured. "
            "Set executor.type in config.yaml or pass --cluster lsf/slurm."
        )

    executor_type = executor_type.lower()
    if executor_type == "lsf":
        return LSFJobManager()
    elif executor_type == "slurm":
        return SlurmJobManager()
    else:
        raise ValueError(
            f"Unknown executor type: {executor_type!r}. "
            "Supported values: 'lsf', 'slurm'."
        )
