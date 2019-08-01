from threading import Lock


class ResourceChangeBuffer:
    def __init__(self):
        self._lock = Lock()
        self._changes = dict()

    def set_changed(self, resource_id, resource_value):
        with self._lock:
            self._changes[resource_id] = resource_value

    def flush_changed(self):
        with self._lock:
            changes = self._changes
            self._changes = dict()
            return changes
