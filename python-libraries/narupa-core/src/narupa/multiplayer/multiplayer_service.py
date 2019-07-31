# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing an implementation of a multiplayer service,.
"""
import logging
import time
from typing import Iterator, Dict

import grpc
import narupa.protocol.multiplayer.multiplayer_pb2 as multiplayer_proto
from narupa.multiplayer import scene
from narupa.multiplayer.multiplayer_lock import MultiplayerObjectLock
from narupa.multiplayer.pubsub import PubSub
from narupa.protocol.multiplayer.multiplayer_pb2 import StreamEndedResponse, Avatar, ResourceRequestResponse, SetResourceValueRequest, CreatePlayerRequest, CreatePlayerResponse, SubscribePlayerAvatarsRequest
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
        default_scene_properties = scene.get_default_scene_properties()
        self.scene_properties = MultiplayerObjectLock(default_scene_properties)
        self.logger = logging.getLogger(__name__)

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
        """
        Subscribe to the avatar stream, infinitely yielding avatars as they are published by other players.
        :param request: Request to join avatar stream.
        :param context: gRPC context.
        :return: yields Avatars
        """
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
        """
        Subscribe to changes in the multiplayer scene properties
        :param request: Request to join scene property changes.
        :param context: gRPC context.
        :return: 
        """
        for publication in self.scene_pubsub.subscribe(request, context):
            yield publication

    def AcquireResourceLock(self,
                            request: multiplayer_proto.AcquireLockRequest,
                            context) -> ResourceRequestResponse:
        """
        Locks the scene for setting properties.
        A scene can only be locked by one user a time, providing their player_id for ownership.
        :param request: Request to lock scene
        :param context: gRPC context.
        :return: Lock scene property indicating whether scene was successfully locked by the client.
        """
        have_lock = self.scene_properties.try_lock(request.player_id)
        self.logger.debug(f'Lock request from player: {request.player_id}. Granted: {have_lock}')
        response = ResourceRequestResponse(success=have_lock)
        return response

    def ReleaseResourceLock(self,
                            request: multiplayer_proto.ReleaseLockRequest,
                            context) -> ResourceRequestResponse:
        """
        Request to unlock the scene.
        A scene can only be unlocked by whoever locked it, until the lock period has elapsed.
        :param request: Request to unlock the scene.
        :param context: gRPC context.
        :return: Reply indicating whether the scene property was unset.
        """
        is_unlocked = self.scene_properties.try_unlock(request.player_id)
        return ResourceRequestResponse(success=is_unlocked)

    def SetResourceValue(self,
                         request: SetResourceValueRequest,
                         context) -> ResourceRequestResponse:
        """
        Sets the scene properties.
        The scene can only be modified if it has been locked by the requester.
        :param request: Request to set scene properties.
        :param context: gRPC context.
        :return: Reply indicating whether scene property set was successful.
        """
        player_id = request.player_id
        resource_id = request.resource_id
        resource_value = request.resource_value

        success = self.scene_properties.try_update_object(player_id, resource_value)
        self.scene_pubsub.publish(request)
        self.logger.debug(f'Scene edit from player: {request.player_id}. Succeeded: {success}')
        return ResourceRequestResponse(success=success)

    def generate_player_id(self):
        """
        Generates a new player ID.

        :return: A unique player ID.
        """
        if len(self.players) == 0:
            return str(1)
        return str(len(self.players) + 1)


