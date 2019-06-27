# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing an implementation of a multiplayer service,.
"""
import logging
from typing import Iterator, Dict

import grpc
import narupa.protocol.multiplayer.multiplayer_pb2 as multiplayer_proto
from narupa.multiplayer import scene
from narupa.multiplayer.multiplayer_lock import MultiplayerObjectLock
from narupa.multiplayer.pubsub import PubSub
from narupa.protocol.multiplayer.multiplayer_pb2 import PublishAvatarReply, Avatar, SceneProperties, EditObjectReply
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import MultiplayerServicer


class MultiplayerService(MultiplayerServicer):
    """
    A generic multiplayer service.
    Provides player IDs for multiplayer sessions, as well as streams for synchronising avatars and the properties
    of the shared multiplayer scene.

    :param max_avatar_queue_size: Maximum number of avatar requests to place in queue.
    :param send_self: Whether to publish updates back to self.

    """
    avatar_pubsub: PubSub
    scene_pubsub: PubSub
    max_queue_size: int
    scene_properties: MultiplayerObjectLock
    players: Dict[object, multiplayer_proto.JoinMultiplayerRequest]

    def __init__(self, max_avatar_queue_size: int = 10, send_self=False):

        super().__init__()

        self.players = {}
        self.avatar_pubsub = PubSub(max_avatar_queue_size)
        self.scene_pubsub = PubSub(max_avatar_queue_size)
        self.max_queue_size = max_avatar_queue_size
        default_scene_properties = scene.get_default_scene_properties()
        self.scene_properties = MultiplayerObjectLock(default_scene_properties)
        self.send_self = send_self
        self.logger = logging.getLogger(__name__)

    def JoinMultiplayer(self, request, context):
        """
        Adds the player to the multiplayer service.
        :param request: Request to join multiplayer.
        :param context: gRPC context.
        :return: A message containing the assigned player ID.
        """
        player_id = self.generate_player_id()
        self.players[player_id] = request
        self.logger.info(f'Username {player_id} has joined multiplayer.')
        return multiplayer_proto.JoinMultiplayerResponse(player_id=player_id)

    def SubscribeToAvatars(self, request, context) -> Avatar:
        """
        Subscribe to the avatar stream, infinitely yielding avatars as they are published by other players.
        :param request: Request to join avatar stream.
        :param context: gRPC context.
        :return: yields Avatars
        """
        if request.player_id not in self.players:
            context.set_code(grpc.StatusCode.INVALID_ARGUMENT)
            context.set_details("Unknown player ID. Join multiplayer session first!")
            self.logger.error('Unknown playerID attempted to join avatars.')
            return
        try:
            for publication in self.avatar_pubsub.subscribe(request, context):
                self.logger.debug(f'Publishing avatar for player: {request.player_id} ')
                yield publication
        except KeyError as e:
            self.logger.error('Unknown playerID used in subscription to avatars.')
            context.set_code(grpc.StatusCode.INVALID_ARGUMENT)
            context.set_details(str(e))
            return

    def PublishAvatar(self, request_iterator: Iterator, context) -> PublishAvatarReply:
        """
        Publish avatar to the avatar stream.
        :param request_iterator: Stream of avatar publications from the client.
        :param context: gRPC context.
        :return: Reply indicating success once stream has ended.
        """

        self.avatar_pubsub.start_publish(request_iterator, context, send_self=self.send_self)
        return PublishAvatarReply()

    def SubscribeToSceneProperties(self, request, context) -> SceneProperties:
        """
        Subscribe to changes in the multiplayer scene properties
        :param request: Request to join scene property changes.
        :param context: gRPC context.
        :return: 
        """
        for publication in self.scene_pubsub.subscribe(request, context):
            yield publication

    def SetLockScene(self, request: multiplayer_proto.LockRequest, context) -> multiplayer_proto.LockRequest:
        """
        Locks the scene for setting properties.
        A scene can only be locked by one user a time, providing their player_id for ownership.
        :param request: Request to lock scene
        :param context: gRPC context.
        :return: Lock scene property indicating whether scene was successfully locked by the client.
        """
        player_id = request.player_id
        is_locked = self.scene_properties.try_lock(player_id)
        reply = multiplayer_proto.LockRequest()
        reply.locked = is_locked
        reply.player_id = player_id
        self.logger.debug(f'Lock request from player: {request.player_id}. Granted: {is_locked}')

        return reply

    def SetSceneProperty(self, request: multiplayer_proto.ScenePropertyRequest, context):
        """
        Sets the scene properties.
        The scene can only be modified if it has been locked by the requester.
        :param request: Request to set scene properties.
        :param context: gRPC context.
        :return: Reply indicating whether scene property set was successful.
        """

        player_id = request.player_id
        properties = request.properties
        success = self.scene_properties.try_update_object(player_id, properties)
        self.scene_pubsub.publish(request, self.send_self)
        reply = EditObjectReply()
        reply.success = success
        reply.player_id = player_id
        self.logger.debug(f'Scene edit from player: {request.player_id}. Succeeded: {success}')

        return reply

    def UnlockScene(self, request, context):
        """
        Request to unlock the scene.
        A scene can only be unlocked by whoever locked it, until the lock period has elapsed.
        :param request: Request to unlock the scene.
        :param context: gRPC context.
        :return: Reply indicating whether the scene property was unset.
        """
        player_id = request.player_id
        is_unlocked = self.scene_properties.try_unlock(player_id)
        reply = multiplayer_proto.LockRequest()
        reply.locked = not is_unlocked
        reply.player_id = player_id
        return reply

    def generate_player_id(self):
        """
        Generates a new player ID.

        :return: A unique player ID.
        """
        if len(self.players) == 0:
            return str(1)
        return str(len(self.players) + 1)


