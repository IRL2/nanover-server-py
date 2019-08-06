from threading import RLock


class ResourceLockedException(Exception):
    pass


class KeyLockableMap:
    def __init__(self):
        self._lock = RLock()
        self._key_lock_owners = dict()
        self._values = dict()

    def player_can_lock_key(self, player_id, key):
        with self._lock:
            return self._key_lock_owners.get(key, player_id) == player_id

    def lock_key(self, owner_id, key):
        with self._lock:
            if not self.player_can_lock_key(owner_id, key):
                raise ResourceLockedException
            self._key_lock_owners[key] = owner_id

    def release_key(self, owner_id, key):
        with self._lock:
            if not self.player_can_lock_key(owner_id, key):
                raise ResourceLockedException
            del self._key_lock_owners[key]

    def remove_owner(self, owner_id):
        with self._lock:
            self._key_lock_owners = {key: owner for key, owner in self._key_lock_owners.items() if owner != owner_id}

    def set(self, owner_id, key, value):
        with self._lock:
            if not self.player_can_lock_key(owner_id, key):
                raise ResourceLockedException
            self._values[key] = value

    def get(self, resource_id, default=None):
        with self._lock:
            return self._values.get(resource_id, default)

    def get_all(self):
        with self._lock:
            return dict(self._values)
