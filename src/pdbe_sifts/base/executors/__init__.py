import os
import subprocess
import time
from abc import ABC, abstractmethod
from collections import defaultdict
from enum import Enum
from textwrap import dedent

from pdbe_sifts.config import load_config
from tqdm import tqdm

from pdbe_sifts.base.log import logger

conf = load_config()


class JobStatus(Enum):
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"
    OTHER = "OTHER"


class Executor(ABC):
    """Executor interface"""

    @abstractmethod
    def submit(
        self,
        job_name,
        command,
        **kwargs,
    ):
        pass

    @abstractmethod
    def get_job_status_by_id(self, job_id) -> str:
        pass

    def monitor(self, frequency):
        """Monitor a list of job IDs."""
        with tqdm(total=len(self.job_ids), desc="Job completion", disable=None) as pbar:
            total_jobs = len(self.job_ids)
            completed_jobs = 0
            cancelled_jobs = 0
            failed_jobs = set()
            while True:
                job_statuses: dict[JobStatus, set] = defaultdict(set)
                for job_id in self.job_ids.copy():
                    status_text = self.get_job_status_by_id(job_id)
                    status = self.map_job_status(status_text)

                    if status == JobStatus.CANCELLED:
                        cancelled_jobs += 1
                        logger.warning(f"Job {job_id} was cancelled. Not requeuing")
                        self.job_ids.remove(job_id)
                    elif status == JobStatus.COMPLETED:
                        completed_jobs += 1
                        self.job_ids.remove(job_id)
                        if job_id in failed_jobs:
                            failed_jobs.remove(job_id)
                    elif status == JobStatus.FAILED:
                        try:
                            logger.debug(f"Job {job_id} failed. Requeuing")
                            self.requeue_job(job_id)
                        except Exception:
                            logger.error(f"Error requeuing job {job_id}", exc_info=True)
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

    @abstractmethod
    def requeue_job(self, job_id):
        pass

    @abstractmethod
    def map_job_status(self):
        pass


class LocalJobManager(Executor):
    procs = []

    def build_command(self, job_name, command, kwargs):
        return command

    def get_job_status_by_id(self, job_id) -> str:
        pass

    def map_job_status(self):
        return super().map_job_status()


class SlurmJobManager(Executor):
    job_ids: set[int] = set()

    def submit(self, job_name, command, **kwargs):
        """Submit a job to Slurm and return the job ID"""
        command = self.build_command(job_name, command, **kwargs)
        logger.debug(command)
        for _ in range(kwargs.get("jobs", 1)):
            output = subprocess.run(command, check=True, capture_output=True)
            job_id = int(output.stdout.decode("utf-8").split()[-1])
            logger.debug(f"Submitted: {job_id}")
            self.job_ids.add(job_id)

        self.monitor(kwargs.get("monitor_frequency", 5))

    def build_command(self, job_name, command, **kwargs):
        log_dir = kwargs.get("log_dir", conf.executor.log_base_dir)
        sbatch_args = [
            "sbatch",
            "--job-name",
            job_name,
            "--partition",
            kwargs.get("queue", conf.executor.default_queue),
            "--cpus-per-task",
            str(kwargs.get("cores", conf.executor.default_cpu_cores)),
            "--gres",
            str(kwargs.get("gpus", "0")),
            "--mem",
            str(kwargs.get("memory", conf.executor.default_memory)),
            "--time",
            str(kwargs.get("timeout_mins") or conf.executor.default_timeout_mins),
            "--chdir",
            log_dir,
            "--open-mode",
            "append",
            "--output",
            kwargs.get("log_file", "slurm-%j.out"),
        ]

        os.makedirs(log_dir, exist_ok=True)
        script_file = os.path.join(log_dir, "slurm.script")
        with open(script_file, "w") as f:
            cmd = f"""\
            #!/bin/bash
            set -e
            {command}
            """
            f.write(dedent(cmd))
            f.flush()

        sbatch_args += [script_file]

        return sbatch_args

    def requeue_job(self, job_id):
        """Requeue a single job."""
        cmd = ["scontrol", "requeue", str(job_id)]
        subprocess.run(cmd, check=True)
        logger.info(f"Requeued failed job {job_id}")

    def get_job_status_by_id(self, job_id) -> str:
        """Get the status of a single job by ID."""
        cmd = [
            "squeue",
            "--noheader",
            "-t",
            "all",
            "--Format",
            "State",
            "--jobs",
            str(job_id),
        ]
        squeue_output = subprocess.run(cmd, check=True, capture_output=True)
        status_text = squeue_output.stdout.decode("utf-8").strip()
        return status_text

    def map_job_status(self, status_text: str) -> JobStatus:
        if status_text in ["FAILED", "TIMEOUT", "OUT_OF_MEMORY"]:
            return JobStatus.FAILED
        elif status_text in ["COMPLETED"]:
            return JobStatus.COMPLETED
        elif status_text in ["RUNNING", "PENDING", "CONFIGURING", "COMPLETING"]:
            return JobStatus.RUNNING
        elif status_text in ["CANCELLED"]:
            return JobStatus.CANCELLED
        else:
            return JobStatus.OTHER


class LSFJobManager(Executor):
    job_ids: set[int] = set()

    def submit(self, job_name, command, **kwargs):
        """Submit a job to LSF and return the job ID"""
        logger.debug(f"Command: {command}")
        logger.info(f"Job name: {job_name}")
        logger.info(f"Job options: {kwargs}")
        cmd = self.build_command(job_name, command, **kwargs)
        logger.debug(f"bsub command: {cmd}")
        for _ in range(kwargs.get("jobs", 1)):
            bsub_output = subprocess.run(cmd, check=True, capture_output=True)
            job_id = int(bsub_output.stdout.decode("utf-8").split()[1][1:-1])
            logger.debug(f"Submitted: {job_id}")
            self.job_ids.add(job_id)

        self.monitor(kwargs.get("monitor_frequency", 5))

    def build_command(self, job_name, command, **kwargs):
        """Builds LSF bsub command based on provided arguments.

        Args:
            job_name (str): Job name
            command (str): Command to execute
            kwargs (dict): Additional arguments
                        Processed arguments are:
                        queue (str): Queue name
                        cpu_cores (int): Number of CPU cores
                        memory (str): Memory limit
                        timeout_mins (int): Timeout in minutes
                        log_dir (str): Log directory
                        log_file (str): Log file
                        error_file (str): Error file

        Returns:
            list: bsub command as list of strings
        """
        bsub_args = [
            "bsub",
            "-J",
            job_name,
            "-q",
            kwargs.get("queue", conf.executor.default_queue),
            "-n",
            str(kwargs.get("cpu_cores", conf.executor.default_cpu_cores)),
            "-M",
            str(kwargs.get("memory", conf.executor.default_memory)),
            "-W",
            str(kwargs.get("timeout_mins") or conf.executor.default_timeout_mins),
            "-cwd",
            kwargs.get("log_dir", conf.executor.log_base_dir),
            "-o",
            kwargs.get("log_file", "lsf-%J.out"),
        ]

        bsub_args += [command]
        return bsub_args

    def requeue_job(self, job_id):
        """Requeue a single job."""
        cmd = ["brequeue", "-e", str(job_id)]
        subprocess.run(cmd, check=True)
        logger.debug(f"Successfully equeued failed job {job_id}")

    def get_job_status_by_id(self, job_id) -> str:
        """Get the status of a single job by ID."""
        cmd = ["bjobs", "-o", "stat", "-noheader", str(job_id)]
        bjobs_output = subprocess.run(cmd, check=True, capture_output=True)
        status_text = bjobs_output.stdout.decode("utf-8").strip()
        if status_text == "EXIT":
            cmd = ["bjobs", "-l", str(job_id)]
            bjobs_output = subprocess.run(cmd, check=True, capture_output=True)
            output = bjobs_output.stdout.decode("utf-8").strip()
            if "TERM_OWNER" in output:
                status_text = "CANCELLED"
        return status_text

    def map_job_status(self, status_text: str) -> JobStatus:
        if status_text in ["DONE"]:
            return JobStatus.COMPLETED
        elif status_text in ["RUN", "PEND", "PSUSP"]:
            return JobStatus.RUNNING
        elif status_text in ["EXIT"]:
            return JobStatus.FAILED
        elif status_text in ["CANCELLED"]:
            return JobStatus.CANCELLED
        else:
            return JobStatus.OTHER


def get_executor(executor=None) -> Executor:
    executor_type = executor or os.getenv("ORC_EXECUTOR") or conf.executor.type
    if not executor_type:
        raise ValueError("No executor specified")
    executor_type = executor_type.lower()
    if executor_type == "slurm":
        return SlurmJobManager()
    elif executor_type == "lsf":
        return LSFJobManager()
    else:
        raise ValueError(
            f"Unknown executor type: {executor_type}. Supported types: slurm, lsf"
        )
