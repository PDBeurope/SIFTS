# Interface for queue managers.
import abc
from typing import Any


class IQueue(abc.ABC):
    queue_type: str

    @abc.abstractmethod
    def push(self, queue_name: str, val: Any):
        pass

    @abc.abstractmethod
    def pop(self, queue_name: str):
        pass

    @abc.abstractmethod
    def pop_and_push(self, source_queue: str, target_queue: str):
        pass

    @abc.abstractmethod
    def length(self, queue_name: str):
        pass

    @abc.abstractmethod
    def get_all(self, queue_name: str):
        pass

    @abc.abstractmethod
    def delete(self, queue_name: str):
        pass

    @abc.abstractmethod
    def rename(self, source_queue: str, target_queue: str):
        pass

    @abc.abstractmethod
    def remove_item(self, queue_name: str, item: Any):
        pass
