"""Singleton factory that creates and caches a queue instance.

Queue parameters are read lazily from the pdbe_sifts OmegaConf config
(``queue.type``, ``queue.host``, ``queue.port``) so that this module can be
imported without needing a config file to be present at import time.

Callers can also pass explicit ``queue_type``, ``host``, and ``port`` to
override config values.
"""

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.queues import IQueue


class BatchableQueue:
    """Singleton factory for :class:`IQueue` instances.

    Only one queue instance is created per (queue_type, host, port) tuple.
    If you call ``BatchableQueue()`` twice with the same parameters you get
    the same underlying queue.

    Usage::

        queue = BatchableQueue().get_queue()
        queue.push("my_queue", "entry_1abc")
    """

    _singleton_queue: IQueue | None = None
    _instance = None
    _config: tuple | None = None

    def __new__(
        cls,
        queue_type: str | None = None,
        host: str | None = None,
        port: int | None = None,
        **kwargs,
    ):
        # Load config lazily (only when a BatchableQueue is actually instantiated).
        from pdbe_sifts.config import load_config

        cfg = load_config()
        queue_type = queue_type or cfg.queue.type
        host = host or cfg.queue.host
        port = int(port or cfg.queue.port)

        # For Redis: read credentials from config if not provided as kwargs
        if queue_type == "redis":
            if "password" not in kwargs:
                cfg_password = getattr(cfg.queue, "password", None)
                if cfg_password:
                    kwargs["password"] = cfg_password
            if "db" not in kwargs:
                cfg_db = getattr(cfg.queue, "db", 0)
                if cfg_db:
                    kwargs["db"] = int(cfg_db)

        new_config = (queue_type, host, port)

        if cls._instance is None or cls._config != new_config:
            logger.info(f"Creating queue manager: type={queue_type} at {host}:{port}")
            cls._instance = super().__new__(cls)
            cls._config = new_config
            cls._singleton_queue = cls._create_queue(queue_type, host, port, **kwargs)
        else:
            logger.debug(f"Re-using existing queue manager at {cls._config}")

        return cls._instance

    def get_queue(self) -> IQueue:
        """Return the singleton queue instance."""
        return self._singleton_queue

    @classmethod
    def _create_queue(
        cls,
        queue_type: str,
        host: str,
        port: int,
        **kwargs,
    ) -> IQueue:
        if queue_type == "local":
            from pdbe_sifts.base.queues.local_queue import LocalQueue

            q = LocalQueue(host=host, port=port)
            q._manager.start()
            return q

        elif queue_type == "redis":
            from pdbe_sifts.base.queues.redis_queue import RedisQueue

            return RedisQueue(
                host=host,
                port=port,
                db=kwargs.get("db", 0),
                password=kwargs.get("password", None),
            )
        else:
            raise ValueError(
                f"Unsupported queue type: {queue_type!r}. "
                "Expected 'local' or 'redis'. Check queue.type in config.yaml."
            )
