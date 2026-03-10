from typing import Any

from redis import StrictRedis
from redis.exceptions import ConnectionError

from pdbe_sifts.base.queues.batchable_queue import IQueue
from pdbe_sifts.config import load_config

conf = load_config()


class RedisQueue(IQueue):
    queue_type = "redis"

    def __init__(self, host: str, port: int, db: int, **kwargs):
        self._server = StrictRedis(
            host=host, port=port, db=db, decode_responses=True, **kwargs
        )
        try:
            self._server.ping()
        except ConnectionError as ce:
            raise Exception("Redis server not running") from ce

    def push(self, queue_name: str, val: Any):
        self._server.rpush(queue_name, val)

    def pop(self, queue_name: str):
        return self._server.rpop(queue_name)

    def pop_and_push(self, source_queue: str, target_queue: str):
        return self._server.rpoplpush(source_queue, target_queue)

    def length(self, queue_name: str):
        return self._server.llen(queue_name)

    def get_all(self, queue_name: str):
        return self._server.lrange(queue_name, 0, -1)

    def delete(self, queue_name: str):
        self._server.delete(queue_name)

    def rename(self, source_queue: str, target_queue: str):
        self._server.rename(source_queue, target_queue)

    def remove_item(self, queue_name: str, item: Any):
        self._server.lrem(queue_name, 1, item)
