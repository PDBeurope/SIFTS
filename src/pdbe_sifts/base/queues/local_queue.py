from collections import defaultdict
from multiprocessing.managers import SyncManager
from typing import Any

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.queues.batchable_queue import IQueue


class QueueManager(SyncManager):
    pass


class LocalQueue(IQueue):
    """A queue that can be used to communicate between processes.

    The queue is implemented using a multiprocessing manager as a list.
    This is to enable remove_item operations on the queue since the default
    multiprocessing queue does not support this.
    """

    _queues: dict[str, list] = defaultdict(list)
    queue_type = "local"

    def __init__(self, host="localhost", port=50000, source_task_id="sifts"):
        logger.info("Created LocalQueue")
        address = (host, port)

        self._manager = QueueManager(address, authkey=source_task_id.encode())
        logger.info(f"Starting queue manager at {address}")

    def push(self, queue_name: str, val: Any):
        """Push a value to the queue. If the queue does not exist, it will be created."""
        self._register(queue_name)
        self._queues[queue_name].append(val)

    def _register(self, queue_name: str):
        if queue_name not in self._queues:
            self._queues[queue_name] = self._manager.list()

    def pop(self, queue_name: str):
        try:
            return self._queues.get(queue_name, []).pop()
        except IndexError:
            return None

    def pop_and_push(self, source_queue: str, target_queue: str):
        try:
            val = self._queues[source_queue].pop()
            self._queues[target_queue].append(val)
            return val
        except IndexError:
            return None

    def length(self, queue_name: str):
        return len(self._queues.get(queue_name, []))

    def get_all(self, queue_name: str):
        return self._queues.get(queue_name, [])

    def delete(self, queue_name: str):
        self._queues.pop(queue_name, None)

    def rename(self, source_queue: str, target_queue: str):
        self._queues[target_queue] = self._queues.pop(source_queue)

    def remove_item(self, queue_name: str, item: Any):
        try:
            self._queues[queue_name].remove(item)
        except ValueError:
            pass
