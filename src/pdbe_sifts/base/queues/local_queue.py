"""Local multiprocessing queue backend.

Uses :class:`multiprocessing.managers.SyncManager` to expose a shared list
over TCP so that both the main process and remote LSF worker processes (which
must reach the login node on the configured port) can share a queue.
"""

from collections import defaultdict
from multiprocessing.managers import SyncManager
from typing import Any

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.queues import IQueue


class QueueManager(SyncManager):
    pass


class LocalQueue(IQueue):
    """Queue backed by a :class:`multiprocessing.managers.SyncManager`.

    The manager exposes a shared list over TCP so that multiple processes
    (including remote LSF/Slurm jobs connecting over the network) can
    enqueue and dequeue entries.

    Args:
        host: Hostname the manager server listens on.
        port: TCP port the manager server listens on.
        source_task_id: Used as the ``authkey`` for the manager (bytes-encoded).
    """

    _queues: dict[str, list] = defaultdict(list)
    queue_type = "local"

    def __init__(
        self,
        host: str = "localhost",
        port: int = 50000,
        source_task_id: str = "pdbe_sifts",
    ):
        logger.info("Created LocalQueue")
        address = (host, port)
        self._manager = QueueManager(address, authkey=source_task_id.encode())
        logger.info(f"Starting queue manager at {address}")

    def push(self, queue_name: str, val: Any) -> None:
        """Push *val* onto *queue_name*, creating the queue if it does not exist."""
        self._register(queue_name)
        self._queues[queue_name].append(val)

    def _register(self, queue_name: str) -> None:
        if queue_name not in self._queues:
            self._queues[queue_name] = self._manager.list()

    def pop(self, queue_name: str) -> Any:
        try:
            return self._queues.get(queue_name, []).pop()
        except IndexError:
            return None

    def pop_and_push(self, source_queue: str, target_queue: str) -> Any:
        try:
            val = self._queues[source_queue].pop()
            self._queues[target_queue].append(val)
            return val
        except IndexError:
            return None

    def length(self, queue_name: str) -> int:
        return len(self._queues.get(queue_name, []))

    def get_all(self, queue_name: str) -> list:
        return list(self._queues.get(queue_name, []))

    def delete(self, queue_name: str) -> None:
        self._queues.pop(queue_name, None)

    def rename(self, source_queue: str, target_queue: str) -> None:
        self._queues[target_queue] = self._queues.pop(source_queue)

    def remove_item(self, queue_name: str, item: Any) -> None:
        try:
            self._queues[queue_name].remove(item)
        except ValueError:
            pass
