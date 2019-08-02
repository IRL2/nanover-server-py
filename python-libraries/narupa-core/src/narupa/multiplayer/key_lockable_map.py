from threading import RLock


class ResourceLockedException(Exception):
    pass


class KeyLockableMap:
    def __init__(self):
        self._lock = RLock()
        self._locks = dict()
        self._values = dict()

    def can_lock(self, player_id, key):
        with self._lock:
            return self._locks.get(key, player_id) == player_id

    def lock_key(self, owner_id, key):
        with self._lock:
            if not self.can_lock(owner_id, key):
                raise ResourceLockedException
            self._locks[key] = owner_id

    def release_key(self, owner_id, key):
        with self._lock:
            if not self.can_lock(owner_id, key):
                raise ResourceLockedException
            del self._locks[key]

    def remove_owner(self, owner_id):
        with self._lock:
            self._locks = {key: owner for key, owner in self._locks.items() if owner != owner_id}

    def set(self, owner_id, key, value):
        with self._lock:
            if not self.can_lock(owner_id, key):
                raise ResourceLockedException
            self._values[key] = value

    def get(self, resource_id, default=None):
        with self._lock:
            return self._values.get(resource_id, default)

    def get_all(self):
        with self._lock:
            return dict(self._values)
