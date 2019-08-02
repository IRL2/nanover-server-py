# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing an implementation of a multiplayer service,.
"""
import logging
from typing import Iterator

import narupa.protocol.multiplayer.multiplayer_pb2 as multiplayer_proto
from narupa.multiplayer.dictionary_change_buffer import DictionaryChangeMultiView, ObjectClosedException
from narupa.multiplayer.key_lockable_map import KeyLockableMap, ResourceLockedException
from narupa.protocol.multiplayer.multiplayer_pb2 import StreamEndedResponse, Avatar, ResourceRequestResponse, SetResourceValueRequest, CreatePlayerRequest, CreatePlayerResponse, SubscribePlayerAvatarsRequest, ResourceValuesUpdate
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import MultiplayerServicer


class MultiplayerService(MultiplayerServicer):
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
        player_id = self.generate_player_id()
        self.players[player_id] = request
        self.logger.info(f'{request.player_name} ({player_id}) has joined multiplayer.')
        return CreatePlayerResponse(player_id=player_id)

    def SubscribePlayerAvatars(self,
                               request: SubscribePlayerAvatarsRequest,
                               context) -> Avatar:
        change_buffer = self._avatars.create_view()
        while context.is_active():
            try:
                changes = change_buffer.flush_changed_blocking()
                for player_id, avatar in changes.items():
                    if player_id != request.ignore_player_id:
                        yield avatar
                #time.sleep(request.update_interval)
            except ObjectClosedException:
                pass

    def UpdatePlayerAvatar(self,
                           request_iterator: Iterator,
                           context) -> StreamEndedResponse:
        for avatar in request_iterator:
            self._avatars.update({avatar.player_id: avatar})
        return StreamEndedResponse()

    def SubscribeAllResourceValues(self, request, context):
        change_buffer = self._resources.create_view()
        while context.is_active():
            try:
                changes = change_buffer.flush_changed_blocking()
                response = ResourceValuesUpdate()
                for key, value in changes.items():
                    entry = response.resource_value_changes.get_or_create(key)
                    entry.MergeFrom(value)
                yield response
                # time.sleep(request.update_interval)  # TODO: give change buffers a wait function
            except ObjectClosedException:
                break

    def AcquireResourceLock(self,
                            request: multiplayer_proto.AcquireLockRequest,
                            context) -> ResourceRequestResponse:
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
        try:
            self.resources.release_key(request.player_id, request.resource_id)
            success = True
        except ResourceLockedException:
            success = False
        return ResourceRequestResponse(success=success)

    def SetResourceValue(self,
                         request: SetResourceValueRequest,
                         context) -> ResourceRequestResponse:
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


