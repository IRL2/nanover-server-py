# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Reference multiplayer client implementation.

"""

from typing import Callable, Sequence

import grpc
import narupa.protocol.multiplayer.multiplayer_pb2_grpc as mult_proto_grpc
from narupa.core import NarupaStubClient
from narupa.multiplayer.multiplayer_server import DEFAULT_PORT
from narupa.protocol.multiplayer.multiplayer_pb2 import (
    CreatePlayerRequest,
    SetResourceValueRequest,
    RemoveResourceKeyRequest,
    AcquireLockRequest,
    ReleaseLockRequest,
    SubscribeAllResourceValuesRequest,
)
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import MultiplayerStub
from narupa.utilities.protobuf_utilities import object_to_value, struct_to_dict

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
        self._value_update_callbacks = set()
        self.resources = dict()

    def close(self):
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
        response = self.stub.CreatePlayer(CreatePlayerRequest(player_name=player_name), timeout=5)
        self._player_id = response.player_id
        if join_streams:
            self.subscribe_all_value_updates()

        return self._player_id

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
        request = SubscribeAllResourceValuesRequest(update_interval=interval)
        self.threads.submit(self._join_scene_properties, request)

    def try_lock_resource(self, resource_id, duration=0):
        """
        Attempt to gain exclusive write access to a particular key in the
        shared key/value store.

        :param resource_id: Key to lock.
        :param duration: Duration to lock the key for (0 for indefinite).
        """
        lock_request = AcquireLockRequest(
            player_id=self.player_id,
            resource_id=resource_id,
            timeout_duration=duration,
        )
        reply = self.stub.AcquireResourceLock(lock_request)
        return reply.success

    def try_release_resource(self, resource_id):
        """
        Attempt to release exclusive write access of a particular key in the
        shared key/value store.

        :param resource_id: Key to release.
        """
        lock_request = ReleaseLockRequest(
            player_id=self.player_id,
            resource_id=resource_id,
        )
        reply = self.stub.ReleaseResourceLock(lock_request)
        return reply.success

    def try_set_resource_value(self, resource_id, value) -> bool:
        """
        Attempt to write a value to a key in the shared key/value store.

        :param resource_id: Key to write to.
        :param value: Value to write.
        """

        value = object_to_value(value)

        request = SetResourceValueRequest(
            player_id=self.player_id,
            resource_id=resource_id,
            resource_value=value,
        )
        response = self.stub.SetResourceValue(request)
        return response.success

    def try_remove_resource_key(self, resource_id: str) -> bool:
        """
        Attempt to remove a key from the shared key/value store.

        :param resource_id: Key to remove.
        """
        request = RemoveResourceKeyRequest(
            player_id=self.player_id,
            resource_id=resource_id,
        )
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
    def _join_scene_properties(self, request):
        for update in self.stub.SubscribeAllResourceValues(request):
            self.resources.update(struct_to_dict(update.resource_value_changes))
            for key in update.resource_value_removals:
                self.resources.pop(key, None)
            keys = set(update.resource_value_changes.keys())
            for callback in self._value_update_callbacks:
                callback(keys)
