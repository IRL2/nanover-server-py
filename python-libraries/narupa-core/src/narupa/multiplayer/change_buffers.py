from typing import Callable, Optional
from contextlib import contextmanager
from threading import Lock, Condition
from narupa.core.timing import yield_interval


class ObjectFrozenException(Exception):
    """
    Raised when an operation on an object cannot be performed because the
    object has been frozen.
    """
    pass


class DictionaryChangeMultiView:
    """
    Provides a means to acquire multiple independent DictionaryChangeBuffers
    tracking a shared dictionary.
    """
    def __init__(self):
        self._content = {}
        self._frozen = False
        self._lock = Lock()
        self._views = set()

    @contextmanager
    def create_view(self):
        """
        Returns a new DictionaryChangeBuffer that tracks changes to the
        shared dictionary, starting with the initial values.
        """
        with self._lock:
            view = DictionaryChangeBuffer()
            view.update(self._content)
            if self._frozen:
                view.freeze()
            self._views.add(view)
        yield view
        self.remove_view(view)

    def remove_view(self, view):
        """
        Freeze the given change buffer and stop providing updates to it.
        """
        with self._lock:
            self._views.remove(view)
            view.freeze()

    def subscribe_changes(self, interval: float = 0):
        """
        Iterates over changes to the shared dictionary, starting with the
        initial values. Waits at least :interval: seconds between each
        iteration.
        """
        with self.create_view() as view:
            yield from view.subscribe_changes(interval)

    def update(self, updates):
        """
        Updates the shared dictionary with key values pairs from :updates:.
        """
        with self._lock:
            if self._frozen:
                raise ObjectFrozenException()
            self._content.update(updates)
            for view in set(self._views):
                try:
                    view.update(updates)
                except ObjectFrozenException:
                    self._views.remove(view)

    def freeze(self):
        """
        Prevent any further updates to the shared dictionary, ensuring that
        future views and subscriptions are frozen and provide a single update
        with the final values.
        """
        with self._lock:
            self._frozen = True
            for view in self._views:
                view.freeze()


class DictionaryChangeBuffer:
    """
    Tracks the latest values of keys that have changed between checks.
    """
    def __init__(self):
        self._frozen = False
        self._lock = Lock()
        self._any_changes = Condition(self._lock)
        self._changes = {}

    def freeze(self):
        """
        Freeze the buffer, ensuring that it cannot be updated with any more
        changes.

        It will still be possible to flush changes one last time if there were
        any unflushed at the time of freezing.
        """
        with self._lock:
            self._frozen = True
            self._any_changes.notify()

    def update(self, updates):
        """
        Update the known changes from a dictionary of keys that have changed
        to their new values.
        """
        with self._lock:
            if self._frozen:
                raise ObjectFrozenException()
            self._changes.update(updates)
            self._any_changes.notify()

    def flush_changed_blocking(self):
        """
        Wait until there are changes and then return them, clearing all
        tracked changes.
        """
        with self._any_changes:
            while not self._changes:
                if self._frozen:
                    raise ObjectFrozenException()
                self._any_changes.wait()
            changes = self._changes
            self._changes = dict()
            return changes

    def subscribe_changes(self, interval: float = 0):
        """
        Iterates over changes to the buffer. Waits at least :interval: seconds
        between each iteration.
        """
        for dt in yield_interval(interval):
            try:
                yield self.flush_changed_blocking()
            except ObjectFrozenException:
                break
