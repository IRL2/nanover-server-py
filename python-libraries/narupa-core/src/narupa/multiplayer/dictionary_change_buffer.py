from threading import Lock, Condition


class ObjectClosedException(Exception):
    pass


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
