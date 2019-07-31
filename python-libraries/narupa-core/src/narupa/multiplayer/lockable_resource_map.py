from threading import RLock


class ResourceLockedException(Exception):
    pass


class LockableResourceMap:
    def __init__(self):
        self._lock = RLock()
        self._locks = dict()
        self._values = dict()

    def lock_key(self, owner_id, key):
        with self._lock:
            if self._locks.get(key, owner_id) != owner_id:
                raise ResourceLockedException
            self._locks[key] = owner_id

    def release_key(self, owner_id, key):
        with self._lock:
            if self._locks.get(key, owner_id) != owner_id:
                raise ResourceLockedException
            del self._locks[key]

    def set(self, owner_id, key, value):
        with self._lock:
            if self._locks.get(key, owner_id) != owner_id:
                raise ResourceLockedException
            self._values[key] = value

    def get(self, resource_id, default=None):
        with self._lock:
            return self._values.get(resource_id, default)
