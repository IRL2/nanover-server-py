# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing an implementation of the :class:`StateServicer`.

"""
from typing import Iterable

from narupa.core.grpc_utils import (
    subscribe_rpc_termination,
    RpcAlreadyTerminatedError,
)
from narupa.core.key_lockable_map import ResourceLockedException

from narupa.multiplayer.change_buffers import (
    DictionaryChange,
)

from narupa.protocol.state import (
    StateServicer,
    StateUpdate,
    SubscribeStateUpdatesRequest,
    UpdateStateRequest,
    UpdateLocksRequest,
    UpdateStateResponse,
    UpdateLocksResponse,
)
from narupa.state.state_dictionary import StateDictionary


class StateService(StateServicer):
    _state_dictionary: StateDictionary

    def __init__(self):
        super().__init__()
        self._id = "service"
        self._state_dictionary = StateDictionary()

    def SubscribeStateUpdates(
            self,
            request: SubscribeStateUpdatesRequest,
            context
    ) -> Iterable[StateUpdate]:
        """
        Provides a stream of updates to a shared key/value store.
        """
        interval = request.update_interval
        with self._state_dictionary.create_view() as change_buffer:
            try:
                subscribe_rpc_termination(context, change_buffer.freeze)
            except RpcAlreadyTerminatedError:
                return
            for change in change_buffer.subscribe_changes(interval):
                yield dictionary_change_to_state_update(change)

    def UpdateState(
            self,
            request: UpdateStateRequest,
            context
    ) -> UpdateStateResponse:
        """
        Attempts an atomic update of the shared key/value store. If any key
        cannot be updates, no change will be made.
        """
        success = True
        change = state_update_to_dictionary_change(request.update)
        try:
            self._state_dictionary.update_state(request.access_token, change)
        except ResourceLockedException:
            success = False
        return UpdateStateResponse(success)

    def UpdateLocks(
            self,
            request: UpdateLocksRequest,
            context,
    ) -> UpdateLocksResponse:
        """
        Attempts to acquire and release locks on keys in the shared key/value
        store. If any of the locks cannot be acquired, none of them will be.
        """
        success = True
        acquire, release = locks_update_to_dictionary_change(request)
        try:
            self._state_dictionary.update_locks(request.access_token, acquire, release)
        except ResourceLockedException:
            success = False
        return UpdateLocksResponse(success)


def state_update_to_dictionary_change(update: StateUpdate) -> DictionaryChange:
    changes = {}
    removals = set()

    for key, value in update:
        if value is None:
            removals.add(key)
        else:
            changes[key] = value

    return DictionaryChange(changes, removals)


def dictionary_change_to_state_update(change: DictionaryChange) -> StateUpdate:
    changes, removals = change

    update = StateUpdate()
    update.changed_keys.update(changes)
    for key in removals:
        update.changed_keys[key] = None

    return update


def locks_update_to_dictionary_change(update: UpdateLocksRequest) -> DictionaryChange:
    release = set()
    acquire = {}

    for key, duration in update.lock_keys:
        if duration > 0:
            acquire[key] = duration
        else:
            release.add(key)

    return DictionaryChange(acquire, release)
