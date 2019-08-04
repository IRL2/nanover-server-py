# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing an implementation of a multiplayer service,.
"""
import logging
from typing import Iterator

import narupa.protocol.multiplayer.multiplayer_pb2 as multiplayer_proto
from narupa.multiplayer.change_buffers import DictionaryChangeMultiView
from narupa.multiplayer.key_lockable_map import KeyLockableMap, ResourceLockedException
from narupa.protocol.multiplayer.multiplayer_pb2 import StreamEndedResponse, Avatar, ResourceRequestResponse, SetResourceValueRequest, CreatePlayerRequest, CreatePlayerResponse, SubscribePlayerAvatarsRequest, ResourceValuesUpdate
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import MultiplayerServicer


class MultiplayerService(MultiplayerServicer):
    """Implementation of the Multiplayer service."""
    def __init__(self):
        super().__init__()

        self.players = {}
        self.logger = logging.getLogger(__name__)

        self._avatars = DictionaryChangeMultiView()
        self._resources = DictionaryChangeMultiView()
        self.resources = KeyLockableMap()

    def CreatePlayer(self,
                     request: CreatePlayerRequest,
                     context) -> CreatePlayerResponse:
        """Create a new unique player and return their id."""
        player_id = self.generate_player_id()
        self.players[player_id] = request
        self.logger.info(f'{request.player_name} ({player_id}) has joined multiplayer.')
        return CreatePlayerResponse(player_id=player_id)

    def SubscribePlayerAvatars(self,
                               request: SubscribePlayerAvatarsRequest,
                               context) -> Avatar:
        """Provides a stream of updates to player avatars."""
        for changes in self._avatars.subscribe_updates(request.update_interval):
            for player_id, avatar in changes.items():
                if player_id != request.ignore_player_id:
                    yield avatar

    def UpdatePlayerAvatar(self,
                           request_iterator: Iterator,
                           context) -> StreamEndedResponse:
        """Accepts a stream of avatar updates."""
        for avatar in request_iterator:
            self._avatars.update({avatar.player_id: avatar})
        return StreamEndedResponse()

    def SubscribeAllResourceValues(self, request, context):
        """Provides a stream of updates to a shared key/value store."""
        for changes in self._resources.subscribe_updates(interval=request.update_interval):
            response = ResourceValuesUpdate()
            for key, value in changes.items():
                entry = response.resource_value_changes.get_or_create(key)
                entry.MergeFrom(value)
            yield response

    def AcquireResourceLock(self,
                            request: multiplayer_proto.AcquireLockRequest,
                            context) -> ResourceRequestResponse:
        """Attempt to acquire exclusive write access to a key in the shared
        key/value store."""
        try:
            self.resources.lock_key(request.player_id, request.resource_id)
            success = True
        except ResourceLockedException:
            success = False
        self.logger.debug(f'{request.player_id} attempts lock on {request.resource_id} (Success: {success})')
        response = ResourceRequestResponse(success=success)
        return response

    def ReleaseResourceLock(self,
                            request: multiplayer_proto.ReleaseLockRequest,
                            context) -> ResourceRequestResponse:
        """Attempt to release exclusive write access to a key in the shared
        key/value store."""
        try:
            self.resources.release_key(request.player_id, request.resource_id)
            success = True
        except ResourceLockedException:
            success = False
        return ResourceRequestResponse(success=success)

    def SetResourceValue(self,
                         request: SetResourceValueRequest,
                         context) -> ResourceRequestResponse:
        """Attempt to write a value in the shared key/value store."""
        try:
            self.resources.set(request.player_id,
                               request.resource_id,
                               request.resource_value)
            success = True
            self._resources.update({request.resource_id: request.resource_value})
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
        self._avatars.close()
        self._resources.close()
