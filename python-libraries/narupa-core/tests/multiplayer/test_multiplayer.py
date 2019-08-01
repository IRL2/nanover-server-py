# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Integration tests of the multiplayer server with the reference multiplayer client.
"""

import time
import pytest

from narupa.multiplayer.multiplayer_client import MultiplayerClient
from narupa.multiplayer.multiplayer_server import MultiplayerServer
from narupa.protocol.multiplayer.multiplayer_pb2 import Avatar, AvatarComponent
import narupa.protocol.multiplayer.multiplayer_pb2 as mult_proto
from google.protobuf.struct_pb2 import Value, Struct


@pytest.fixture
def multiplayer_server():
    server = MultiplayerServer(address='localhost', port=0)
    yield server
    server.close()


@pytest.fixture
def multiplayer_server_send_self():
    server = MultiplayerServer(address='localhost', port=0)
    yield server
    server.close()


@pytest.fixture
def multiplayer_server_client(multiplayer_server):
    client = MultiplayerClient(port=multiplayer_server.port)
    yield multiplayer_server, client
    client.close()


@pytest.fixture
def avatar():
    components = [AvatarComponent(name="Head",
                                  position=[0, 0, 1],
                                  rotation=[1, 1, 1, 1])]
    avatar = Avatar(player_id="1", component=components)
    return avatar


def test_join_multiplayer(multiplayer_server_client):
    multiplayer_server, multiplayer_client = multiplayer_server_client
    player_id = multiplayer_client.join_multiplayer("user", join_streams=False)
    print('Player ID: ', player_id)
    assert player_id == "1"


def test_join_multiplayer_twice(multiplayer_server_client):
    multiplayer_server, multiplayer_client = multiplayer_server_client
    multiplayer_client.join_multiplayer("user", join_streams=False)
    player_id = multiplayer_client.join_multiplayer("user")
    assert player_id == "1"


def test_join_avatar_stream(multiplayer_server_client):
    multiplayer_server, multiplayer_client = multiplayer_server_client
    multiplayer_client.join_avatar_stream()


def test_join_avatar_before_joining_multiplayer(multiplayer_server_client):
    multiplayer_server, multiplayer_client = multiplayer_server_client
    request = mult_proto.SubscribePlayerAvatarsRequest()
    avatar_stream = multiplayer_client.stub.SubscribePlayerAvatars(request)


def test_join_avatar_stream_twice(multiplayer_server_client):
    multiplayer_server, multiplayer_client = multiplayer_server_client
    multiplayer_client.join_multiplayer(player_name="user", join_streams=False)

    request = mult_proto.SubscribePlayerAvatarsRequest()
    avatar_stream1 = multiplayer_client.stub.SubscribePlayerAvatars(request)
    avatar_stream2 = multiplayer_client.stub.SubscribePlayerAvatars(request)


def test_join_publish_avatar(multiplayer_server_client):
    multiplayer_server, multiplayer_client = multiplayer_server_client
    multiplayer_client.join_multiplayer(player_name="user")
    multiplayer_client.join_avatar_publish()


def test_publish_avatar_to_self(multiplayer_server_client, avatar):
    multiplayer_server, multiplayer_client = multiplayer_server_client
    player_id = multiplayer_client.join_multiplayer(player_name="user")
    time.sleep(0.05)

    multiplayer_client.join_avatar_publish()
    multiplayer_client.join_avatar_stream(ignore_self=False)
    time.sleep(0.05)

    avatar.player_id = player_id
    multiplayer_client.publish_avatar(avatar)
    time.sleep(0.05)

    assert str(multiplayer_client.current_avatars[player_id]) == str(avatar)


def test_publish_avatar_ignore_self(multiplayer_server_client, avatar):
    multiplayer_server, multiplayer_client = multiplayer_server_client
    player_id = multiplayer_client.join_multiplayer(player_name="user")
    time.sleep(0.05)

    multiplayer_client.join_avatar_publish()
    multiplayer_client.join_avatar_stream(ignore_self=True)
    time.sleep(0.05)

    avatar.player_id = player_id
    multiplayer_client.publish_avatar(avatar)
    time.sleep(0.05)

    assert len(multiplayer_client.current_avatars) == 0


def test_publish_avatar_multiple_transmission(multiplayer_server_client, avatar):
    multiplayer_server, multiplayer_client = multiplayer_server_client
    player_id = multiplayer_client.join_multiplayer(player_name="user")
    time.sleep(0.05)

    multiplayer_client.join_avatar_publish()
    multiplayer_client.join_avatar_stream(ignore_self=False)
    time.sleep(0.05)

    multiplayer_client.publish_avatar(avatar)
    avatar.component[0].position[:] = [0, 0, 2]
    multiplayer_client.publish_avatar(avatar)
    avatar.component[0].position[:] = [0, 0, 3]
    multiplayer_client.publish_avatar(avatar)
    time.sleep(0.05)

    client_avatar = multiplayer_client.current_avatars[player_id]
    assert client_avatar.component[0].position == [0, 0, 3]


@pytest.fixture
def test_scene():
    pose = Struct()
    pose["position"] = {"x": 1, "y": 1, "z": 1}
    pose["rotation"] = {"x": 0, "y": 0, "z": 0, "w": 1}
    pose["scale"] = 1
    return Value(struct_value=pose)


def test_can_lock_unlocked(multiplayer_server_client):
    server, client = multiplayer_server_client
    assert client.try_lock_resource("scene")


def test_can_lock_own_locked(multiplayer_server_client):
    server, client = multiplayer_server_client
    client.try_lock_resource("scene")
    assert client.try_lock_resource("scene")


def test_can_release_own_lock(multiplayer_server_client):
    server, client = multiplayer_server_client
    client.try_lock_resource("scene")
    assert client.try_release_resource("scene")


def test_can_set_unlocked(multiplayer_server_client, test_scene):
    server, client = multiplayer_server_client
    assert client.try_set_resource_value("scene", test_scene)


def test_can_set_own_locked(multiplayer_server_client, test_scene):
    server, client = multiplayer_server_client
    client.try_lock_resource("scene")
    assert client.try_set_resource_value("scene", test_scene)


def test_set_value_updates_server_values(multiplayer_server_client, test_scene):
    server, client = multiplayer_server_client
    client.try_set_resource_value("scene", test_scene)
    server_scene = server.multiplayer_service.resources.get("scene")
    assert str(test_scene) == str(server_scene)


def test_set_value_sends_update(multiplayer_server_client, test_scene):
    server, client = multiplayer_server_client
    client.subscribe_all_value_updates()
    time.sleep(0.1)
    client.try_set_resource_value("scene", test_scene)
    time.sleep(0.1)
    recv_scene = client.resources.get("scene")
    assert str(test_scene) == str(recv_scene)


def test_server_sends_initial_values(multiplayer_server_client, test_scene):
    server, client = multiplayer_server_client
    client.try_set_resource_value("scene", test_scene)
    time.sleep(0.1)
    assert client.resources.get("scene") is None
    client.subscribe_all_value_updates()
    time.sleep(0.1)
    assert str(client.resources.get("scene")) == str(test_scene)


def test_cant_lock_other_locked(multiplayer_server_client):
    server, client = multiplayer_server_client
    client.try_lock_resource("scene", player_id="fake")
    assert not client.try_lock_resource("scene")


def test_cant_release_other_lock(multiplayer_server_client):
    server, client = multiplayer_server_client
    client.try_lock_resource("scene", player_id="fake")
    assert not client.try_release_resource("scene")


def test_cant_set_other_locked(multiplayer_server_client, test_scene):
    server, client = multiplayer_server_client
    client.try_lock_resource("scene", player_id="fake")
    assert not client.try_set_resource_value("scene", test_scene)
