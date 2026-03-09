"""Queue interface and factory for cluster job coordination.

Queues are used in cluster (LSF/Slurm) mode to distribute entries across
compute nodes.  Two backends are provided:

  * LocalQueue  — multiprocessing.managers.SyncManager over TCP (default)
  * RedisQueue  — Redis-backed queue (requires ``redis`` extra)

The factory :class:`BatchableQueue` reads the queue type from the
pdbe_sifts config (``queue.type``) and returns a singleton instance.
"""

import abc
from typing import Any


class IQueue(abc.ABC):
    """Abstract interface for queue managers used by cluster jobs."""

    queue_type: str

    @abc.abstractmethod
    def push(self, queue_name: str, val: Any) -> None:
        """Append *val* to the tail of *queue_name*."""

    @abc.abstractmethod
    def pop(self, queue_name: str) -> Any:
        """Remove and return the tail element of *queue_name*, or ``None``."""

    @abc.abstractmethod
    def pop_and_push(self, source_queue: str, target_queue: str) -> Any:
        """Atomically pop from *source_queue* and push onto *target_queue*.

        Returns the moved element, or ``None`` if *source_queue* is empty.
        """

    @abc.abstractmethod
    def length(self, queue_name: str) -> int:
        """Return the number of elements in *queue_name*."""

    @abc.abstractmethod
    def get_all(self, queue_name: str) -> list:
        """Return all elements of *queue_name* without removing them."""

    @abc.abstractmethod
    def delete(self, queue_name: str) -> None:
        """Delete *queue_name* entirely."""

    @abc.abstractmethod
    def rename(self, source_queue: str, target_queue: str) -> None:
        """Rename *source_queue* to *target_queue*."""

    @abc.abstractmethod
    def remove_item(self, queue_name: str, item: Any) -> None:
        """Remove the first occurrence of *item* from *queue_name*."""
