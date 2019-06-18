"""
Provides a dictionary of queues.
"""

from typing import Dict, Hashable, Generator
from queue import Queue
from threading import Lock
from contextlib import contextmanager


class DictOfQueues:
    """
    Dictionary of request queues.

    .. example::

        # A thread working on its own queue
        many_queues = DictOfQueues()
        request_id = 0
        with many_queues.one_queue(request_id) as queue:
            # do stuff

        # A thread populating all the queues
        many_queues = DictOfQueues(queue_max_size=100)
        for queue in many_queue.iter_queue():
            queue.put('something')


    Adding or removing a key from the dictionary *must* be done while holding
    :attr:`DictOfQueues.lock`. The recommended way to register and unregister
    a queue is to use the :meth:`one_queue` context manager.
    """
    queue_max_size: int
    queues: Dict[Hashable, Queue]
    lock: Lock

    def __init__(self, queue_max_size=0):
        """
        :param queue_max_size: The maximum size of each queue. If set to 0,
            there is no maximum size and queues can grow indefinitely. When a
            queue reaches its maximum size, adding a element is blocking until
            the size of the queue decreases.
        """
        self.queue_max_size = queue_max_size
        self.queues = {}
        self.lock = Lock()

    @contextmanager
    def one_queue(self, request_id):
        """
        Works with a queue.

        This method is a context manager that creates and registers the queue,
        provides the queue to the calling scope, and un-registers the queue
        when exiting the context.

        :param request_id: The key for the queue. This key has to be unique, if
            a queue is already registered with that key, then a
            :exc:`ValueError` is raised.
        """
        queue = Queue(maxsize=self.queue_max_size)
        with self.lock:
            if request_id in self.queues:
                raise ValueError(f'The key {request_id} is already registered.')
            self.queues[request_id] = queue
        try:
            yield queue
        finally:
            with self.lock:
                del self.queues[request_id]

    def iter_queues(self) -> Generator[Queue, None, None]:
        """
        Iterate over the queues.

        The method places a lock on the dictionary so no queue can be added or
        removed while iterating.
        """
        with self.lock:
            yield from self.queues.values()
