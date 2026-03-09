"""Redis queue backend.

Requires the ``redis`` extra::

    pip install pdbe_sifts[redis]
    # or: pip install redis
"""

from typing import Any

from pdbe_sifts.base.queues import IQueue


class RedisQueue(IQueue):
    """Queue backed by a Redis server.

    Args:
        host: Redis server hostname.
        port: Redis server port.
        db: Redis database number (default 0).
        kwargs: Additional keyword arguments forwarded to
            :class:`redis.StrictRedis`.
    """

    queue_type = "redis"

    def __init__(self, host: str, port: int, db: int = 0, **kwargs):
        try:
            from redis import StrictRedis
            from redis.exceptions import ConnectionError as RedisConnectionError
        except ImportError as exc:
            raise ImportError(
                "The 'redis' package is required for RedisQueue. "
                "Install it with: pip install pdbe_sifts[redis]"
            ) from exc

        self._server = StrictRedis(
            host=host, port=port, db=db, decode_responses=True, **kwargs
        )
        try:
            self._server.ping()
        except RedisConnectionError as ce:
            raise RuntimeError(
                f"Could not connect to Redis at {host}:{port} — "
                "is the Redis server running?"
            ) from ce

    def push(self, queue_name: str, val: Any) -> None:
        self._server.rpush(queue_name, val)

    def pop(self, queue_name: str) -> Any:
        return self._server.rpop(queue_name)

    def pop_and_push(self, source_queue: str, target_queue: str) -> Any:
        return self._server.rpoplpush(source_queue, target_queue)

    def length(self, queue_name: str) -> int:
        return self._server.llen(queue_name)

    def get_all(self, queue_name: str) -> list:
        return self._server.lrange(queue_name, 0, -1)

    def delete(self, queue_name: str) -> None:
        self._server.delete(queue_name)

    def rename(self, source_queue: str, target_queue: str) -> None:
        self._server.rename(source_queue, target_queue)

    def remove_item(self, queue_name: str, item: Any) -> None:
        self._server.lrem(queue_name, 1, item)
