import time
from threading import Lock, Condition


def yield_interval(interval):
    """
    Yield at a set interval, accounting for the time spent outside of this
    function.
    :param interval: Number of seconds to put between yields
    """
    last_yield = time.monotonic() - interval
    while True:
        time_since_yield = time.monotonic() - last_yield
        wait_duration = max(0, interval - time_since_yield)
        time.sleep(wait_duration)
        yield time.monotonic() - last_yield
        last_yield = time.monotonic()


class ObjectClosedException(Exception):
    """Raised when an operation on an object cannot be performed because the
    object has been closed."""
    pass


class DictionaryChangeMultiView:
    """Provides a means to acquire multiple independent DictionaryChangeBuffers
    tracking a shared dictionary."""
    def __init__(self):
        self._content = {}
        self._closed = False
        self._lock = Lock()
        self._views = set()

    def create_view(self):
        """Returns a new DictionaryChangeBuffer that tracks changes to the
        shared dictionary, starting with the initial values."""
        with self._lock:
            view = DictionaryChangeBuffer()
            view.update(self._content)
            if self._closed:
                view.close()
            self._views.add(view)
            return view

    def subscribe_updates(self, interval=0):
        """Iterates over changes to the shared dictionary, starting with the
        initial values. Waits at least :interval: seconds between each
        iteration."""
        view = self.create_view()
        for dt in yield_interval(interval):
            try:
                yield view.flush_changed_blocking()
            except ObjectClosedException:
                break

    def update(self, updates):
        """Updates the shared dictionary with key values pairs from :updates:.
        """
        with self._lock:
            if self._closed:
                raise ObjectClosedException()
            self._content.update(updates)
            for view in self._views:
                view.update(updates)

    def close(self):
        """Prevent any further updates to the shared dictionary, ensuring that
        future views and subscriptions provide a single update with the final
        values."""
        with self._lock:
            self._closed = True
            for view in self._views:
                view.close()


class DictionaryChangeBuffer:
    """Tracks the latest values of keys that have changed between checks."""
    def __init__(self):
        self._closed = False
        self._lock = Lock()
        self._any_changes = Condition(self._lock)
        self._changes = {}

    def close(self):
        """Close the buffer, ensuring that no more changes will be tracked."""
        with self._lock:
            self._closed = True
            self._any_changes.notify()

    def update(self, updates):
        """Update the known changes from a dictionary of keys that have changed
        to their new values."""
        with self._lock:
            if self._closed:
                raise ObjectClosedException()
            self._changes.update(updates)
            self._any_changes.notify()

    def flush_changed_blocking(self):
        """Wait until there are changes and then return them, clearing all
        tracked changes."""
        with self._any_changes:
            while not self._changes:
                if self._closed:
                    raise ObjectClosedException()
                self._any_changes.wait()
            changes = self._changes
            self._changes = dict()
            return changes
