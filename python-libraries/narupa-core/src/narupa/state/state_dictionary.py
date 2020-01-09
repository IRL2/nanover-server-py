# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from contextlib import contextmanager
from threading import Lock
from typing import ContextManager, Set, Dict

from narupa.core.key_lockable_map import KeyLockableMap, ResourceLockedException
from narupa.multiplayer.change_buffers import (
    DictionaryChangeMultiView,
    DictionaryChangeBuffer,
    DictionaryChange,
)


class StateDictionary:
    _lock: Lock
    _change_views: DictionaryChangeMultiView
    _write_locks: KeyLockableMap

    def __init__(self):
        self._lock = Lock()
        self._change_views = DictionaryChangeMultiView()
        self._write_locks = KeyLockableMap()

    @property
    def content(self) -> Dict[str, object]:
        with self._lock:
            return self._change_views.copy_content()

    def get_change_buffer(self) -> ContextManager[DictionaryChangeBuffer]:
        """
        Return a DictionaryChangeBuffer that tracks changes to this dictionary.
        """
        return self._change_views.create_view()

    def update_state(
            self,
            access_token: object,
            change: DictionaryChange,
    ):
        """
        Update the dictionary with key changes and removals, using any locks
        permitted by the given access token.
        """
        with self._lock:
            if not self._can_token_access_keys(access_token, change.updates.keys()):
                raise ResourceLockedException

            self._change_views.update(change.updates, change.removals)

    def update_locks(
            self,
            access_token: object,
            acquire: Dict[str, float],
            release: Set[str],
    ):
        """
        Acquire and release locks for the given access token.
        """
        with self._lock:
            if not self._can_token_access_keys(access_token, acquire.keys()):
                raise ResourceLockedException

            for key in release:
                try:
                    self._write_locks.release_key(access_token, key)
                except ResourceLockedException:
                    pass # don't care if we can't release a lock
            for key, duration in acquire.items():
                self._write_locks.lock_key(access_token, key, duration)

    def _can_token_access_keys(self, access_token, keys):
        """
        Return whether or not all keys are either unlocked or locked by the
        given access token.
        """
        return all(self._write_locks.player_can_lock_key(access_token, key) for key in keys)
