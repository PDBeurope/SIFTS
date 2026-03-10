from pdbe_sifts.base.log import logger
from pdbe_sifts.base.queues.iqueue import IQueue
from pdbe_sifts.config import load_config

conf = load_config()


class BatchableQueue:
    """Factory class for creating batchable queues."""

    _singleton_queue = None
    _instance = None
    _config = (conf.queue.type, conf.queue.host, conf.queue.port)

    def __new__(
        cls,
        queue_type: str = conf.queue.type,
        host: str = conf.queue.host,
        port: int = conf.queue.port,
        **kwargs,
    ):
        """Initialize a queue based on the config file.

        Returns:
            IQueue: A queue object.
        """

        # Create a singleton queue manager
        _config = (queue_type, host, port)
        cls._queue_type = queue_type
        if cls._instance is None or cls._config != _config:
            logger.info(f"Creating multiprocessing queue manager at {_config}")

            cls._instance = super().__new__(cls)
            cls._config = _config
            cls._singleton_queue = cls._create_queue(**kwargs)
        else:
            logger.debug(f"Using existing queue manager at {cls._config}")

        return cls._instance

    def get_queue(self):
        return self._singleton_queue

    @classmethod
    def _create_queue(cls, **kwargs) -> IQueue:
        queue_type, host, port = cls._config

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
            raise Exception(
                "Unsupported queue type. Check conf.queue.type in config.yaml"
            )
