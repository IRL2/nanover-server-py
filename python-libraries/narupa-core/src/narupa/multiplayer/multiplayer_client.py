# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Reference multiplayer client implementation.

"""

from typing import Dict, Callable, Sequence, Optional

import grpc

import narupa.protocol.multiplayer.multiplayer_pb2 as mult_proto
import narupa.protocol.multiplayer.multiplayer_pb2_grpc as mult_proto_grpc
from narupa.core import NarupaStubClient, NarupaClient
from narupa.utilities.protobuf_utilities import object_to_value, struct_to_dict
from narupa.utilities.request_queues import SingleItemQueue
from narupa.utilities.change_buffers import yield_interval
from narupa.multiplayer.multiplayer_server import DEFAULT_PORT
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import MultiplayerStub

UpdateCallback = Callable[[Sequence[str]], None]


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


class MultiplayerClient(NarupaStubClient):
    """
    Represents a client to the multiplayer server.

    """
    stub: mult_proto_grpc.MultiplayerStub
    current_avatars: Dict[int, mult_proto.Avatar]
    _avatar_queue: SingleItemQueue
    DEFAULT_CONNECTION_PORT = DEFAULT_PORT

    def __init__(self, *,
                 channel: grpc.Channel,
                 make_channel_owner: bool = False,
                 send_interval: float = 1 / 60):
        super().__init__(channel=channel,
                         stub=MultiplayerStub,
                         make_channel_owner=make_channel_owner)
        self._player_id = None
        self._send_interval = send_interval
        self._avatar_queue = SingleItemQueue()
        self._value_update_callbacks = set()
        self.current_avatars = {}
        self.resources = dict()

    def close(self):
        self._avatar_queue.put(None)
        super().close()

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
            self.subscribe_all_value_updates()

        return self._player_id

    def join_avatar_stream(self, interval=0, ignore_self=True):
        """
        Joins the avatar stream, which will start receiving avatar updates in
        the background.

        :param interval: Minimum time (in seconds) between receiving two updates
        for the same player's avatar.
        :param ignore_self: Whether to request the server not to send our own
        avatar updates back to us.
        """
        ignore = self.player_id if ignore_self else None
        request = mult_proto.SubscribePlayerAvatarsRequest(ignore_player_id=ignore,
                                                           update_interval=interval)
        self.threads.submit(self._join_avatar_stream, request)

    def join_avatar_publish(self):
        """
        Joins the avatar publishing stream.

        Use :func:`~MultiplayerClient.publish_avatar` to publish.
        """
        self._assert_has_player_id()
        self.threads.submit(self._join_avatar_publish)

    def publish_avatar(self, avatar: mult_proto.Avatar):
        """
        Updates an avatar to be published.

        :param avatar: The avatar to be published.
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

    def subscribe_all_value_updates(self, interval=0):
        """
        Begin receiving updates to the shared key/value store.

        :param interval: Minimum time (in seconds) between receiving new updates
            for any and all values.
        """
        request = mult_proto.SubscribeAllResourceValuesRequest(update_interval=interval)
        self.threads.submit(self._join_scene_properties, request)

    def try_lock_resource(self, resource_id, duration=0):
        """
        Attempt to gain exclusive write access to a particular key in the
        shared key/value store.

        :param resource_id: Key to lock.
        :param duration: Duration to lock the key for (0 for indefinite).
        """
        lock_request = mult_proto.AcquireLockRequest(player_id=self.player_id,
                                                     resource_id=resource_id,
                                                     timeout_duration=duration)
        reply = self.stub.AcquireResourceLock(lock_request)
        return reply.success

    def try_release_resource(self, resource_id):
        """
        Attempt to release exclusive write access of a particular key in the
        shared key/value store.

        :param resource_id: Key to release.
        """
        lock_request = mult_proto.ReleaseLockRequest(player_id=self.player_id,
                                                     resource_id=resource_id)
        reply = self.stub.ReleaseResourceLock(lock_request)
        return reply.success

    def try_set_resource_value(self, resource_id, value) -> bool:
        """
        Attempt to write a value to a key in the shared key/value store.

        :param resource_id: Key to write to.
        :param value: Value to write.
        """

        value = object_to_value(value)

        request = mult_proto.SetResourceValueRequest(player_id=self.player_id,
                                                     resource_id=resource_id,
                                                     resource_value=value)
        response = self.stub.SetResourceValue(request)
        return response.success

    def try_remove_resource_key(self, resource_id: str) -> bool:
        """
        Attempt to remove a key from the shared key/value store.

        :param resource_id: Key to remove.
        """
        request = mult_proto.RemoveResourceKeyRequest(player_id=self.player_id,
                                                      resource_id=resource_id)
        response = self.stub.RemoveResourceKey(request)
        return response.success

    def add_value_update_callback(self, callback: UpdateCallback):
        """
        Add a callback method to be called whenever a value changes in the
        shared key/value store.

        :param callback: Method to be called, that takes the set of keys that
            have changed as arguments.
        """
        self._value_update_callbacks.add(callback)

    def _assert_has_player_id(self):
        if not self.joined_multiplayer:
            raise RuntimeError("Join multiplayer before attempting this operation.")

    @_end_upon_channel_close
    def _join_avatar_stream(self, request):
        for avatar in self.stub.SubscribePlayerAvatars(request):
            self.current_avatars[avatar.player_id] = avatar

    @_end_upon_channel_close
    def _join_avatar_publish(self):
        response = self.stub.UpdatePlayerAvatar(self._publish_avatar())

    def _publish_avatar(self):
        for dt in yield_interval(self._send_interval):
            avatar = self._avatar_queue.get(block=True)
            if avatar is None:
                break
            yield avatar

    @_end_upon_channel_close
    def _join_scene_properties(self, request):
        for update in self.stub.SubscribeAllResourceValues(request):
            self.resources.update(struct_to_dict(update.resource_value_changes))
            for key in update.resource_value_removals:
                self.resources.pop(key, None)
            keys = set(update.resource_value_changes.keys())
            for callback in self._value_update_callbacks:
                callback(keys)
