# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing an implementation of a multiplayer service,.
"""
import logging
from threading import Lock
from typing import Iterator

import narupa.protocol.multiplayer.multiplayer_pb2 as multiplayer_proto
from narupa.multiplayer.change_buffers import DictionaryChangeMultiView
from narupa.multiplayer.key_lockable_map import KeyLockableMap, ResourceLockedException
from narupa.protocol.multiplayer.multiplayer_pb2 import StreamEndedResponse, Avatar, ResourceRequestResponse, SetResourceValueRequest, CreatePlayerRequest, CreatePlayerResponse, SubscribePlayerAvatarsRequest, ResourceValuesUpdate
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import MultiplayerServicer


class MultiplayerService(MultiplayerServicer):
    """
    Implementation of the Multiplayer service.
    """
    def __init__(self):
        super().__init__()

        self.players = {}
        self.logger = logging.getLogger(__name__)

        self._avatars = DictionaryChangeMultiView()
        self._resource_write_lock = Lock()
        self._resources = DictionaryChangeMultiView()
        self.resources = KeyLockableMap()

    def CreatePlayer(self,
                     request: CreatePlayerRequest,
                     context) -> CreatePlayerResponse:
        """
        Create a new unique player and return their id.
        """
        player_id = self.generate_player_id()
        self.players[player_id] = request
        self.logger.info(f'{request.player_name} ({player_id}) has joined '
                         f'multiplayer.')
        return CreatePlayerResponse(player_id=player_id)

    def SubscribePlayerAvatars(self,
                               request: SubscribePlayerAvatarsRequest,
                               context) -> Avatar:
        """
        Provides a stream of updates to player avatars.
        """
        interval = request.update_interval
        with self._avatars.create_view() as change_buffer:
            if not context.add_callback(lambda: change_buffer.freeze()):
                return
            for changes in change_buffer.subscribe_changes(interval):
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
        with self._resources.create_view() as change_buffer:
            if not context.add_callback(lambda: change_buffer.freeze()):
                return
            for changes in change_buffer.subscribe_changes(interval):
                response = ResourceValuesUpdate()
                for key, value in changes.items():
                    entry = response.resource_value_changes.get_or_create(key)
                    entry.MergeFrom(value)
                yield response

    def AcquireResourceLock(self,
                            request: multiplayer_proto.AcquireLockRequest,
                            context) -> ResourceRequestResponse:
        """
        Attempt to acquire exclusive write access to a key in the shared
        key/value store.
        """
        try:
            duration = request.timeout_duration
            if duration <= 0:
                duration = None
            self.resources.lock_key(request.player_id,
                                    request.resource_id,
                                    duration)
            success = True
        except ResourceLockedException:
            success = False
        self.logger.debug(f'{request.player_id} attempts lock on {request.resource_id} (Success: {success})')
        response = ResourceRequestResponse(success=success)
        return response

    def ReleaseResourceLock(self,
                            request: multiplayer_proto.ReleaseLockRequest,
                            context) -> ResourceRequestResponse:
        """
        Attempt to release exclusive write access to a key in the shared
        key/value store.
        """
        try:
            self.resources.release_key(request.player_id, request.resource_id)
            success = True
        except ResourceLockedException:
            success = False
        return ResourceRequestResponse(success=success)

    def SetResourceValue(self,
                         request: SetResourceValueRequest,
                         context) -> ResourceRequestResponse:
        """
        Attempt to write a value in the shared key/value store.
        """
        try:
            # TODO: single lockable+subscribable structure?
            with self._resource_write_lock:
                self.resources.set(request.player_id,
                                   request.resource_id,
                                   request.resource_value)
                self._resources.update({request.resource_id: request.resource_value})
            success = True
        except ResourceLockedException:
            success = False

        self.logger.debug(f'{request.player_id} attempts {request.resource_id}={request.resource_value} (Successs: {success})')
        return ResourceRequestResponse(success=success)

    def generate_player_id(self):
        """
        Generates a new player ID.

        :return: A unique player ID.
        """
        return str(len(self.players) + 1)

    def close(self):
        self._avatars.freeze()
        self._resources.freeze()

    def _clear_player_avatar(self, player_id: str):
        """
        Replace a player's avatar with an empty avatar.
        """
        avatar = Avatar(player_id=player_id, component=[])
        self._avatars.update({avatar.player_id: avatar})
