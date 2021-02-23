# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing utility classes used by the multiplayer service to create a
shared key/value store between multiple clients.
"""
from contextlib import contextmanager
from threading import Lock, Condition
from typing import Any, Set, Dict, Iterator, Iterable, Optional, Generator, Mapping

from .timing import yield_interval

KeyUpdates = Mapping[str, Any]
KeyRemovals = Iterable[str]


class DictionaryChange:
    updates: KeyUpdates
    removals: KeyRemovals

    def __init__(
            self,
            updates: Optional[KeyUpdates] = None,
            removals: Optional[KeyRemovals] = None,
    ):
        self.updates = updates or {}
        self.removals = removals or set()

    def __iter__(self):
        return iter((self.updates, self.removals))


class ObjectFrozenError(Exception):
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
    _content: Dict[str, Any]
    _frozen: bool
    _lock: Lock
    _views: Set['DictionaryChangeBuffer']

    def __init__(self):
        self._content = {}
        self._frozen = False
        self._lock = Lock()
        self._views = set()

    @contextmanager
    def create_view(self) -> Generator['DictionaryChangeBuffer', None, None]:
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

    def remove_view(self, view: 'DictionaryChangeBuffer'):
        """
        Freeze the given change buffer and stop providing updates to it.
        """
        with self._lock:
            self._views.remove(view)
            view.freeze()

    def subscribe_changes(self, interval: float = 0) \
            -> Iterator[DictionaryChange]:
        """
        Iterates over changes to the shared dictionary, starting with the
        initial values. Waits at least :interval: seconds between each
        iteration.
        """
        with self.create_view() as view:
            yield from view.subscribe_changes(interval)

    def update(self, updates: KeyUpdates = None, removals: KeyRemovals = None):
        """
        Updates the shared dictionary with key values pairs from :updates: and
        key removals from :removals:.
        """
        updates = updates or {}
        removals = removals or set()
        with self._lock:
            if self._frozen:
                raise ObjectFrozenError()
            self._update_content(updates, removals)
            self._update_views(updates, removals)

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

    def copy_content(self) -> Dict[str, Any]:
        """
        Return a shallow copy of the dictionary at this instant.
        """
        with self._lock:
            return dict(self._content)

    def _update_content(self, updates: KeyUpdates, removals: KeyRemovals):
        """
        Apply updates and removals to the known content of this dictionary.
        """
        self._content.update(updates)
        for key in removals:
            self._content.pop(key, None)

    def _update_views(self, updates: KeyUpdates, removals: KeyRemovals):
        """
        Forward updates and removals to the tracked change buffers of this
        dictionary.
        """
        for view in set(self._views):
            try:
                view.update(updates, removals)
            except ObjectFrozenError:
                self._views.remove(view)


class DictionaryChangeBuffer:
    """
    Tracks the latest values of keys that have changed between checks.
    """
    _frozen: bool
    _lock: Lock
    _any_changes: Condition
    _pending_changes: Dict[str, Any]
    _pending_removals: Set[str]

    def __init__(self):
        self._frozen = False
        self._lock = Lock()
        self._any_changes = Condition(self._lock)
        self._pending_changes = {}
        self._pending_removals = set()

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

    def update(self, updates: KeyUpdates = None, removals: KeyRemovals = None):
        """
        Update the known changes from a dictionary of keys that have changed
        to their new values or have been removed.
        """
        updates = updates or {}
        removals = removals or set()
        with self._lock:
            if self._frozen:
                raise ObjectFrozenError()

            # pending key deletions with new values are no longer removed
            for key in updates:
                self._pending_removals.discard(key)
            # pending changes for deleted keys don't matter
            for key in removals:
                self._pending_changes.pop(key, None)
            self._pending_changes.update(updates)
            self._pending_removals.update(removals)
            self._any_changes.notify()

    def flush_changed_blocking(self) -> DictionaryChange:
        """
        Wait until there are changes and then return them, clearing all
        tracked changes.
        """
        with self._any_changes:
            while not self._pending_changes and not self._pending_removals:
                if self._frozen:
                    raise ObjectFrozenError()
                self._any_changes.wait()
            changes, removals = self._pending_changes, self._pending_removals
            self._pending_changes, self._pending_removals = dict(), set()
            return DictionaryChange(changes, removals)

    def subscribe_changes(self, interval: float = 0) \
            -> Iterator[DictionaryChange]:
        """
        Iterates over changes to the buffer. Waits at least :interval: seconds
        between each iteration.
        """
        for dt in yield_interval(interval):
            try:
                yield self.flush_changed_blocking()
            except ObjectFrozenError:
                break
