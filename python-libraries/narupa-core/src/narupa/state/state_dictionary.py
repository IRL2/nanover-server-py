# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing `StateDictionary`, a class for tracking and making changes to a
shared key/value store.
"""
from contextlib import contextmanager
from threading import Lock
from typing import ContextManager, Generator, Dict, Optional, Iterable

from narupa.utilities.change_buffers import (
    DictionaryChangeMultiView,
    DictionaryChangeBuffer,
    DictionaryChange,
)
from narupa.utilities.event import Event
from narupa.utilities.key_lockable_map import (
    KeyLockableMap, ResourceLockedError,
)
from narupa.utilities.protobuf_utilities import Serializable


class StateDictionary:
    """
    Mechanism for tracking and making changes to a shared key/value store,
    including the facility to acquire exclusive write access to values.
    """
    content_updated: Event

    _lock: Lock
    _change_views: DictionaryChangeMultiView
    _write_locks: KeyLockableMap

    def __init__(self):
        self.content_updated = Event()
        self._lock = Lock()
        self._change_views = DictionaryChangeMultiView()
        self._write_locks = KeyLockableMap()

    def freeze(self):
        self._change_views.freeze()

    def copy_content(self):
        """
        Return a shallow copy of the dictionary content at this instant.
        """
        with self.lock_content() as content:
            return dict(content)

    @contextmanager
    def lock_content(self) -> Generator[Dict[str, Serializable], None, None]:
        """
        Context manager for reading the current state while delaying any changes
        to it via an exclusive lock.
        """
        with self._lock:
            yield self._change_views.copy_content()

    def get_change_buffer(self) -> ContextManager[DictionaryChangeBuffer]:
        """
        Return a DictionaryChangeBuffer that tracks changes to this dictionary.
        """
        return self._change_views.create_view()

    def update_state(
            self,
            access_token: Serializable,
            change: DictionaryChange,
    ):
        """
        Update the dictionary with key changes and removals, using any locks
        permitted by the given access token. If any key cannot be updated, no
        change will be made.

        :raises ResourceLockedError: if the access token cannot acquire all keys
            for updating.
        """
        with self._lock:
            keys = set(change.updates.keys()) | set(change.removals)
            if not self._can_token_access_keys(access_token, keys):
                raise ResourceLockedError
            self._change_views.update(change.updates, change.removals)
        self.content_updated.invoke(access_token=access_token, change=change)

    def update_locks(
            self,
            access_token: Serializable,
            acquire: Optional[Dict[str, Optional[float]]] = None,
            release: Optional[Iterable[str]] = None,
    ):
        """
        Acquire and release locks for the given access token. If any of the
        locks cannot be acquired, none of them will be. Requested lock releases
        are carried out regardless.

        :raises ResourceLockedError: if the access token cannot acquire all
            requested keys.
        """
        acquire = acquire or {}
        release = release or ()
        with self._lock:
            if not self._can_token_access_keys(access_token, acquire.keys()):
                raise ResourceLockedError

            for key in release:
                try:
                    self._write_locks.release_key(access_token, key)
                except ResourceLockedError:
                    pass  # don't care if we can't release a lock
            for key, duration in acquire.items():
                self._write_locks.lock_key(access_token, key, duration)

    def _can_token_access_keys(self, access_token, keys):
        """
        Return whether or not all keys are either unlocked or locked by the
        given access token.
        """
        return all(self._write_locks.player_can_lock_key(access_token, key) for key in keys)