# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Provides a dictionary of queues.
"""

from typing import Dict, Hashable, Generator, Tuple
from queue import Queue, Empty
from threading import Lock, Condition
from contextlib import contextmanager
from time import monotonic as time

from narupa.protocol.trajectory import GetFrameResponse, FrameData


class DictOfQueues:
    """
    Dictionary of request queues.
    This class is used by Narupa servers to provide a thread-safe
    way to publish data to multiple clients using queues.

    .. code-block:: python

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

    def __contains__(self, key):
        with self.lock:
            return key in self.queues

    @contextmanager
    def one_queue(self, request_id, queue_class=Queue):
        """
        Works with a queue.

        This method is a context manager that creates and registers the queue,
        provides the queue to the calling scope, and un-registers the queue
        when exiting the context.

        :param request_id: The key for the queue. This key has to be unique, if
            a queue is already registered with that key, then a
            :exc:`ValueError` is raised.
        :param queue_class: The class to instantiate for that queue. By default,
            a :class:`Queue` is instantiated.
        """
        queue = queue_class(maxsize=self.queue_max_size)
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

    def iter_queues_items(self) -> Generator[Tuple[Hashable, Queue], None, None]:
        """
        Iterate over the queues and their keys.

        The method places a lock on the dictionary so no queue can be added or
        removed while iterating.
        """
        with self.lock:
            yield from self.queues.items()


# adapted from https://github.com/python/cpython/blob/master/Lib/queue.py
class SingleItemQueue:
    """
    Mimics the basic interface of a :class:`Queue` but only stores one item.
    """

    def __init__(self, maxsize=None):
        """
        :param maxsize: Unused parameter, included for compatibility with
            :class:`Queue`.
        """
        self._lock = Lock()
        self._item = None
        self._has_item = False

        self.not_empty = Condition(self._lock)

    def put(self, item, **kwargs):
        """
        Store a value, replace the previous one if any.

        This method is thread-safe and is meant to be a drop in replacement
        to :meth:`Queue.put`.

        :param item: The value to store.
        :param kwargs: Unused arguments for compatibility with :meth:`Queue.put`.
        """
        with self._lock:
            self._item = item
            self._has_item = True
            self.not_empty.notify()

    def get(self, block=True, timeout=None):
        """
        Get the stored value, and remove it from storage.

        If there is no value to get, then the method raises an :exc:`Empty`
        exception.

        This method is thread-safe and is meant to be a drop in replacement
        to :meth:`Queue.get`.

        :param block: Whether to wait until a value is available.
        :param timeout: Timeout for waiting until a value is available.
        :return: The stored value.
        """
        with self.not_empty:
            if not block:
                if not self._has_item:
                    raise Empty
            elif timeout is None:
                while not self._has_item:
                    self.not_empty.wait()
            elif timeout < 0:
                raise ValueError("'timeout' must be a non-negative number")
            else:
                endtime = time() + timeout
                while not self._has_item:
                    remaining = endtime - time()
                    if remaining <= 0.0:
                        raise Empty
                    self.not_empty.wait(remaining)
            item = self._item
            self._item = None
            self._has_item = False
            return item


class GetFrameResponseAggregatingQueue(SingleItemQueue):
    """
    SingleItemQueue specifically for GetFrameResponse items. Put items will be
    aggregated with any existing item so that there is at most one item in the
    queue at any time.
    """
    def put(self, item: GetFrameResponse, **kwargs):
        with self._lock:
            if item is None:
                # None is the sentinel value to indicate that the queue user
                # should terminate, so it is safe to discard aggregated frames.
                self._item = None
            else:
                if self._item is None:
                    self._item = GetFrameResponse(
                        frame_index=item.frame_index,
                        frame=FrameData(),
                    )
                self._item.frame_index = item.frame_index
                self._item.frame.MergeFrom(item.frame)
            self._has_item = True
            self.not_empty.notify()

