# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing an implementation of a multiplayer service,.
"""
import logging
import time
from threading import Lock
from typing import Iterator

import narupa.protocol.multiplayer.multiplayer_pb2 as multiplayer_proto
from narupa.multiplayer.resource_change_buffer import ResourceChangeBuffer
from narupa.multiplayer.lockable_resource_map import LockableResourceMap, ResourceLockedException
from narupa.multiplayer.pubsub import PubSub
from narupa.protocol.multiplayer.multiplayer_pb2 import StreamEndedResponse, Avatar, ResourceRequestResponse, SetResourceValueRequest, CreatePlayerRequest, CreatePlayerResponse, SubscribePlayerAvatarsRequest, ResourceValuesUpdate
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import MultiplayerServicer


class MultiplayerService(MultiplayerServicer):
    """
    A generic multiplayer service.
    Provides player IDs for multiplayer sessions, as well as streams for synchronising avatars and the properties
    of the shared multiplayer scene.

    :param max_avatar_queue_size: Maximum number of avatar requests to place in queue.
    :param send_self: Whether to publish updates back to self.

    """
    def __init__(self, max_avatar_queue_size: int = 10):
        super().__init__()

        self.players = {}
        self.avatar_pubsub = PubSub(max_avatar_queue_size)
        self.scene_pubsub = PubSub(max_avatar_queue_size)
        self.max_queue_size = max_avatar_queue_size
        self.logger = logging.getLogger(__name__)

        self._change_buffers_lock = Lock()
        self._change_buffers = set()
        self.resources = LockableResourceMap()

    def CreatePlayer(self,
                     request: CreatePlayerRequest,
                     context) -> CreatePlayerResponse:
        """
        Adds the player to the multiplayer service.
        :param request: Request to join multiplayer.
        :param context: gRPC context.
        :return: A message containing the assigned player ID.
        """
        player_id = self.generate_player_id()
        self.players[player_id] = request
        self.logger.info(f'{request.player_name} ({player_id}) has joined multiplayer.')
        return CreatePlayerResponse(player_id=player_id)

    def SubscribePlayerAvatars(self,
                               request: SubscribePlayerAvatarsRequest,
                               context) -> Avatar:
        outgoing_id = str(time.monotonic())
        for publication in self.avatar_pubsub.subscribe(outgoing_id, context):
            self.logger.debug(f'Publishing avatar to stream {outgoing_id}')
            yield publication

    def UpdatePlayerAvatar(self,
                           request_iterator: Iterator,
                           context) -> StreamEndedResponse:
        """
        Publish avatar to the avatar stream.
        :param request_iterator: Stream of avatar publications from the client.
        :param context: gRPC context.
        :return: Reply indicating success once stream has ended.
        """
        self.avatar_pubsub.start_publish(request_iterator, context)
        return StreamEndedResponse()

    def SubscribeAllResourceValues(self, request, context):
        change_buffer = self._create_change_buffer()
        initial_values = ResourceValuesUpdate()
        for key, value in self.resources.get_all().items():
            entry = initial_values.resource_value_changes.get_or_create(key)
            entry.MergeFrom(value)
        yield initial_values

        while True:
            response = ResourceValuesUpdate()
            changes = change_buffer.flush_changed()
            for key, value in changes.items():
                entry = response.resource_value_changes.get_or_create(key)
                entry.MergeFrom(value)
            yield response
            time.sleep(1/30) # TODO: give change buffers a wait function

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
            self._buffer_change(request.resource_id, request.resource_value)
        except ResourceLockedException:
            success = False

        self.logger.debug(f'{request.player_id} attempts {request.resource_id}={request.resource_value} (Successs: {success})')
        return ResourceRequestResponse(success=success)

    def _create_change_buffer(self):
        buffer = ResourceChangeBuffer()
        with self._change_buffers_lock:
            self._change_buffers.add(buffer)
        return buffer

    def _buffer_change(self, key, value):
        with self._change_buffers_lock:
            for buffer in self._change_buffers:
                buffer.set_changed(key, value)

    def generate_player_id(self):
        """
        Generates a new player ID.

        :return: A unique player ID.
        """
        return str(len(self.players) + 1)


