from queue import Empty
from threading import Lock, Condition
from time import monotonic as time

from nanover.trajectory import FrameData


# adapted from https://github.com/python/cpython/blob/master/Lib/queue.py
class LastItemQueue:
    """
    Mimics the basic interface of a :class:`Queue` but only stores one item.
    """

    def __init__(self):
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


class FrameMergingQueue(LastItemQueue):
    """
    SingleItemQueue specifically for GetFrameResponse items. Put items will be
    aggregated with any existing item so that there is at most one item in the
    queue at any time.
    """

    def close(self):
        self.put(None)

    def put(self, item: FrameData | None, **kwargs):
        with self._lock:
            if item is None:
                # None is the sentinel value to indicate that the queue user
                # should terminate, so it is safe to discard aggregated frames.
                self._item = None
            else:
                if self._item is None:
                    self._item = FrameData()
                self._item.update(item)
            self._has_item = True
            self.not_empty.notify()
