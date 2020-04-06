# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing an implementation of a multiplayer service,.
"""
import logging
from typing import Iterator, Callable, Dict, Set, ContextManager

import narupa.protocol.multiplayer.multiplayer_pb2 as multiplayer_proto
from narupa.state.state_dictionary import StateDictionary
from narupa.state.state_service import validate_dict_is_serializable
from narupa.utilities.grpc_utilities import (
    subscribe_rpc_termination,
    RpcAlreadyTerminatedError,
)
from narupa.utilities.protobuf_utilities import value_to_object, \
    deep_copy_serializable_dict
from narupa.utilities.change_buffers import (
    DictionaryChangeMultiView,
    DictionaryChange,
)
from narupa.utilities.key_lockable_map import (
    ResourceLockedError,
)
from narupa.protocol.multiplayer.multiplayer_pb2 import (
    StreamEndedResponse, Avatar, ResourceRequestResponse,
    SetResourceValueRequest, CreatePlayerRequest, CreatePlayerResponse,
    SubscribePlayerAvatarsRequest, ResourceValuesUpdate,
)
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import MultiplayerServicer, add_MultiplayerServicer_to_server

MULTIPLAYER_SERVICE_NAME = "multiplayer"


class MultiplayerService(MultiplayerServicer):
    """
    Implementation of the Multiplayer service.
    """
    _avatars: DictionaryChangeMultiView
    _state_dictionary: StateDictionary

    def __init__(self):
        super().__init__()
        self.name: str = MULTIPLAYER_SERVICE_NAME
        self.add_to_server_method: Callable = add_MultiplayerServicer_to_server

        self.players = {}
        self.logger = logging.getLogger(__name__)

        self._avatars = DictionaryChangeMultiView()
        self._state_dictionary = StateDictionary()

    def lock_state(self) -> ContextManager[Dict[str, object]]:
        """
        Context manager for reading the current state while delaying any changes
        to it.
        """
        return self._state_dictionary.lock_content()

    def copy_state(self) -> Dict[str, object]:
        """
        Return a deep copy of the current state.
        """
        with self.lock_state() as state:
            return deep_copy_serializable_dict(state)

    def update_state(self, access_token: object, change: DictionaryChange):
        """
        Attempts an atomic update of the shared key/value store. If any key
        cannot be updated, no change will be made.

        :raises ResourceLockedError: if the access token cannot acquire all keys
            for updating.
        :raises TypeError: if the update values cannot be serialized for
            transmission.
        """
        validate_dict_is_serializable(change.updates)
        self._state_dictionary.update_state(access_token, change)

    def update_locks(
            self,
            access_token: object,
            acquire: Dict[str, float],
            release: Set[str],
    ):
        """
        Attempts to acquire and release locks on keys in the shared key/value
        store. If any of the locks cannot be acquired, none of them will be.
        Requested lock releases are carried out regardless.

        :raises ResourceLockedError: if the access token cannot acquire all
            requested keys.
        """
        self._state_dictionary.update_locks(access_token, acquire, release)

    def CreatePlayer(self,
                     request: CreatePlayerRequest,
                     context) -> CreatePlayerResponse:
        """
        Create a new unique player and return their id.
        """
        player_id = self.generate_player_id()
        self.players[player_id] = request
        return CreatePlayerResponse(player_id=player_id)

    def SubscribePlayerAvatars(self,
                               request: SubscribePlayerAvatarsRequest,
                               context) -> Avatar:
        """
        Provides a stream of updates to player avatars.
        """
        interval = request.update_interval
        with self._avatars.create_view() as change_buffer:
            try:
                subscribe_rpc_termination(context, change_buffer.freeze)
            except RpcAlreadyTerminatedError:
                return
            for changes, removals in change_buffer.subscribe_changes(interval):
                for player_id, avatar in changes.items():
                    if player_id != request.ignore_player_id:
                        yield avatar

    def UpdatePlayerAvatar(self,
                           request_iterator: Iterator,
                           context) -> StreamEndedResponse:
        """
        Accepts a stream of avatar updates.
        """
        touched_player_ids = set()

        def clear_touched_avatars():
            for player_id in touched_player_ids:
                self._clear_player_avatar(player_id)

        context.add_callback(clear_touched_avatars)

        for avatar in request_iterator:
            touched_player_ids.add(avatar.player_id)
            self._avatars.update({avatar.player_id: avatar})
        return StreamEndedResponse()

    def SubscribeAllResourceValues(self, request, context):
        """
        Provides a stream of updates to a shared key/value store.
        """
        interval = request.update_interval
        with self._state_dictionary.get_change_buffer() as change_buffer:
            try:
                subscribe_rpc_termination(context, change_buffer.freeze)
            except RpcAlreadyTerminatedError:
                return
            for changes, removals in change_buffer.subscribe_changes(interval):
                response = ResourceValuesUpdate()
                response.resource_value_removals.extend(removals)
                response.resource_value_changes.update(changes)
                yield response

    def AcquireResourceLock(self,
                            request: multiplayer_proto.AcquireLockRequest,
                            context) -> ResourceRequestResponse:
        """
        Attempt to acquire exclusive write access to a key in the shared
        key/value store.
        """
        success = True
        try:
            duration = request.timeout_duration
            if duration <= 0:
                duration = None
            self.update_locks(
                request.player_id,
                {request.resource_id: duration},
                set(),
            )
        except ResourceLockedError:
            success = False
        return ResourceRequestResponse(success=success)

    def ReleaseResourceLock(self,
                            request: multiplayer_proto.ReleaseLockRequest,
                            context) -> ResourceRequestResponse:
        """
        Attempt to release exclusive write access to a key in the shared
        key/value store.
        """
        success = True
        try:
            self.update_locks(
                request.player_id,
                {},
                set([request.resource_id]),
            )
        except ResourceLockedError:
            success = False
        return ResourceRequestResponse(success=success)

    def SetResourceValue(
            self,
            request: SetResourceValueRequest,
            context,
    ) -> ResourceRequestResponse:
        """
        Attempt to write a value in the shared key/value store.
        """
        success = True
        resource_value = value_to_object(request.resource_value)
        try:
            self.update_state(
                request.player_id,
                DictionaryChange({request.resource_id: resource_value}, []),
            )
        except ResourceLockedError:
            success = False
        return ResourceRequestResponse(success=success)

    def RemoveResourceKey(
            self,
            request: multiplayer_proto.RemoveResourceKeyRequest,
            context,
    ) -> ResourceRequestResponse:
        """
        Attempt to remove a key from the shared key/value store.
        """
        success = True
        try:
            self.update_state(
                request.player_id,
                DictionaryChange({}, [request.resource_id]),
            )
        except ResourceLockedError:
            success = False
        return ResourceRequestResponse(success=success)

    def generate_player_id(self):
        """
        Generates a new player ID.

        :return: A unique player ID.
        """
        return str(len(self.players) + 1)

    def close(self):
        self._avatars.freeze()

    def _clear_player_avatar(self, player_id: str):
        """
        Replace a player's avatar with an empty avatar.
        """
        avatar = Avatar(player_id=player_id, components=[])
        self._avatars.update({avatar.player_id: avatar})
