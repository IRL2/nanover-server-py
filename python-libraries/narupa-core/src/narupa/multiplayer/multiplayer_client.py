# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Reference multiplayer client implementation.

"""

import time
from queue import Queue, Empty
from typing import Dict
from concurrent import futures

import grpc

from narupa.core import DEFAULT_CONNECT_ADDRESS
from narupa.multiplayer.multiplayer_server import DEFAULT_PORT
import narupa.protocol.multiplayer.multiplayer_pb2 as mult_proto
import narupa.protocol.multiplayer.multiplayer_pb2_grpc as mult_proto_grpc


def _end_upon_channel_close(function):
    """
    Wrapper function to streaming RPCs that gracefully handles any exceptions thrown when the channel has closed.
    :param function:
    :return:
    """
    def wrapped(self, *args, **kwargs):
        try:
            function(self, *args, **kwargs)
        except grpc.RpcError as e:
            if self.closed:
                return
            if e.args and 'Channel closed' in e.args[0]:
                return
            else:
                raise e

    return wrapped


class MultiplayerClient(object):
    """
    Represents a client to the multiplayer server.

    :param address: IP or web address of server to connect to.
    :param port: Server port.
    :param pubsub_fps: FPS at which to send/receive pub/sub messages.
    """
    stub: mult_proto_grpc.MultiplayerStub
    current_avatars: Dict[int, mult_proto.Avatar]
    _avatar_queue: Queue

    def __init__(self, address=DEFAULT_CONNECT_ADDRESS, port=DEFAULT_PORT, pubsub_fps: float = 30, channel=None):

        if channel is None:
            self.channel = grpc.insecure_channel("{0}:{1}".format(address, port))
            self.shared_channel = False
        else:
            self.channel = channel
            self.shared_channel = True
        self.stub = mult_proto_grpc.MultiplayerStub(self.channel)
        self._player_id = None
        self._avatar_queue = Queue()
        self.pubsub_fps = pubsub_fps
        self.pubsub_rate = 1.0 / pubsub_fps
        self.current_avatars = {}
        self.closed = False
        self.scene_properties = None
        self.threadpool = futures.ThreadPoolExecutor(max_workers=10)

    def close(self):
        """
        Closes the client connection.
        :return:
        """
        self.closed = True
        self.threadpool.shutdown(wait=False)
        if not self.shared_channel:
            self.channel.close()

    def join_multiplayer(self, player_name, join_streams=True):
        """
        Joins a multiplayer server
        :param player_name: The user name to use for multiplayer.
        :param join_streams: Whether to automatically join all streams.
        :return: Player ID received from the multiplayer server.
        """
        if self.joined_multiplayer:
            return self._player_id
        response = self.stub.CreatePlayer(mult_proto.CreatePlayerRequest(player_name=player_name), timeout=5)
        self._player_id = response.player_id
        if join_streams:
            self.join_avatar_stream()
            self.join_avatar_publish()
            self.join_scene_properties_stream()

        return self._player_id

    def join_avatar_stream(self):
        """
        Joins the avatar stream, which will start receiving avatar updates in the background.
        """
        self._ensure_joined_multiplayer()
        request = mult_proto.SubscribePlayerAvatarsRequest(ignore_player_id=self.player_id)
        self.threadpool.submit(self._join_avatar_stream, request)

    def join_avatar_publish(self):
        """
        Joins the avatar publishing stream.
        Use publish_avatar to publish.
        """
        self._ensure_joined_multiplayer()
        self.threadpool.submit(self._join_avatar_publish)

    def publish_avatar(self, avatar):
        """
        Updates an avatar to be published.
        :param avatar:
        """
        self._avatar_queue.put(avatar)

    @property
    def player_id(self):
        """
        The player ID assigned to this client after joining multiplayer.
        :return: The player ID.
        """
        return self._player_id

    @property
    def joined_multiplayer(self) -> bool:
        """
        Indicates whether multiplayer joined, and the client has a valid player ID.
        :return: True if multiplayer has been joined, false otherwise.
        """
        return self.player_id is not None

    def join_scene_properties_stream(self):
        self._ensure_joined_multiplayer()
        request = mult_proto.SubscribeAllResourceValuesRequest()
        self.threadpool.submit(self._join_scene_properties_stream, request)

    def try_lock_scene(self):
        lock_request = mult_proto.AcquireLockRequest(player_id=self.player_id,
                                                     resource_id="scene")
        reply = self.stub.AcquireResourceLock(lock_request)
        return reply.success

    def try_unlock_scene(self):
        lock_request = mult_proto.ReleaseLockRequest(player_id=self.player_id,
                                                     resource_id="scene")
        reply = self.stub.ReleaseResourceLock(lock_request)
        return reply.success

    def set_scene_properties(self, properties) -> bool:
        """
        Attempts to set the multiplayer scene properties.
        If the scene properties are locked by someone else, this will fail.
        :param properties: The new properties to apply.
        :return: True if properties sucessfully set, false otherwise.
        """
        self._ensure_joined_multiplayer()
        if not self.try_lock_scene():
            return False
        property_request = mult_proto.SetResourceValueRequest(player_id=self.player_id,
                                                              resource_id="scene",
                                                              resource_value=properties)
        reply = self.stub.SetResourceValue(property_request)
        if not reply.success:
            # This shouldn't happen, as the box will have been locked above, but check it just in case.
            return False  # pragma: no cover
        if not self.try_unlock_scene():
            # Should always be able to unlock the box, but throw exception here just in case.
            raise RuntimeError("Unable to unlock scene after edit!")  # pragma: no cover
        return True

    def _ensure_joined_multiplayer(self):
        if not self.joined_multiplayer:
            raise RuntimeError("Join multiplayer before attempting this operation.")

    @_end_upon_channel_close
    def _join_avatar_stream(self, request):
        for publication in self.stub.SubscribePlayerAvatars(request):
            self.current_avatars[publication.player_id] = publication
            time.sleep(self.pubsub_rate / min(1, len(self.current_avatars)))

    @_end_upon_channel_close
    def _join_avatar_publish(self):
        self.stub.UpdatePlayerAvatar(self._publish_avatar())

    def _publish_avatar(self):
        while True:
            time.sleep(self.pubsub_rate)
            try:
                avatar = self._avatar_queue.get(block=True, timeout=0.5)
            except Empty:
                pass
            else:
                yield avatar

    @_end_upon_channel_close
    def _join_scene_properties(self, request):
        for publication in self.stub.SubscribeAllResourceValues(request):
            self.scene_properties = publication
            time.sleep(self.pubsub_rate)

    @_end_upon_channel_close
    def _join_scene_properties_stream(self, request):
        self.stub.SubscribeAllResourceValues(self._join_scene_properties(request))
