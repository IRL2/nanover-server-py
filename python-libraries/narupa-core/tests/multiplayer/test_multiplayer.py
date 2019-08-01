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
def server_client_pair():
    """
    Provides a multiplayer server hosting on an available port on localhost,
    and a multiplayer client connected to it.
    """
    server = MultiplayerServer(address='localhost', port=0)
    client = MultiplayerClient(port=server.port)

    with client, server:
        yield server, client


@pytest.fixture
def avatar():
    """
    Provides avatar test data.
    """
    components = [AvatarComponent(name="Head",
                                  position=[0, 0, 1],
                                  rotation=[1, 1, 1, 1])]
    avatar = Avatar(player_id="1", component=components)
    return avatar


@pytest.fixture
def scene():
    """
    Provides scene test data.
    """
    pose = Struct()
    pose["position"] = {"x": 1, "y": 1, "z": 1}
    pose["rotation"] = {"x": 0, "y": 0, "z": 0, "w": 1}
    pose["scale"] = 1
    return Value(struct_value=pose)


def test_join_multiplayer(server_client_pair):
    """
    Test that it's possible to join multiplayer and receive a player id.
    """
    server, client = server_client_pair
    player_id = client.join_multiplayer("user", join_streams=False)
    assert player_id is not None


def test_join_multiplayer_twice_same_id(server_client_pair):
    """
    Test that joining multiplayer again gives you your existing player id.
    """
    server, client = server_client_pair
    first_id = client.join_multiplayer("user", join_streams=False)
    second_id = client.join_multiplayer("user", join_streams=False)
    assert first_id == second_id


def test_join_avatar_stream(server_client_pair):
    """
    Test that the avatar stream can be joined.
    """
    server, client = server_client_pair
    client.join_avatar_stream()


def test_cant_publish_avatar_without_player(server_client_pair):
    """
    Test that attempting to join avatar publishing without a player id fails.
    """
    server, client = server_client_pair
    with pytest.raises(RuntimeError):
        client.join_avatar_publish()


def test_join_publish_avatar_with_player(server_client_pair):
    """
    Test that join avatar publish with a player id works.
    """
    server, client = server_client_pair
    client.join_multiplayer(player_name="user")
    client.join_avatar_publish()


def test_publish_avatar_to_self(server_client_pair, avatar):
    """
    Test that a published avatar makes it back to yourself.
    """
    server, client = server_client_pair
    player_id = client.join_multiplayer(player_name="user")
    client.join_avatar_publish()
    client.join_avatar_stream(ignore_self=False)
    time.sleep(0.05)

    avatar.player_id = player_id
    client.publish_avatar(avatar)
    time.sleep(0.05)

    assert str(client.current_avatars[player_id]) == str(avatar)


def test_publish_avatar_ignore_self(server_client_pair, avatar):
    """
    Test that a published avatar does not make it back to yourself if you are
    ignoring yourself.
    """
    server, client = server_client_pair
    player_id = client.join_multiplayer(player_name="user")
    client.join_avatar_publish()
    client.join_avatar_stream()
    time.sleep(0.05)

    avatar.player_id = player_id
    client.publish_avatar(avatar)
    time.sleep(0.05)

    assert len(client.current_avatars) == 0


def test_publish_avatar_multiple_transmission(server_client_pair, avatar):
    """
    Test that multiple different avatar publishes results in the client
    reflecting the last published avatar.
    """
    server, client = server_client_pair
    player_id = client.join_multiplayer(player_name="user")
    client.join_avatar_publish()
    client.join_avatar_stream(ignore_self=False)
    time.sleep(0.05)

    client.publish_avatar(avatar)
    avatar.component[0].position[:] = [0, 0, 2]
    client.publish_avatar(avatar)
    avatar.component[0].position[:] = [0, 0, 3]
    client.publish_avatar(avatar)
    time.sleep(0.05)

    client_avatar = client.current_avatars[player_id]
    assert client_avatar.component[0].position == [0, 0, 3]


def test_can_lock_unlocked(server_client_pair):
    """
    Test that an unlocked resource can be locked.
    """
    server, client = server_client_pair
    assert client.try_lock_resource("scene")


def test_can_lock_own_locked(server_client_pair):
    """
    Test that an attempt to lock a resource you have already lock succeeds.
    """
    server, client = server_client_pair
    client.try_lock_resource("scene")
    assert client.try_lock_resource("scene")


def test_can_release_own_lock(server_client_pair):
    """
    Test that you can release a resource you locked.
    """
    server, client = server_client_pair
    client.try_lock_resource("scene")
    assert client.try_release_resource("scene")


def test_can_set_unlocked(server_client_pair, scene):
    """
    Test that you can set a resource value if it is unlocked.
    """
    server, client = server_client_pair
    assert client.try_set_resource_value("scene", scene)


def test_set_unlocked_repeated(server_client_pair, scene):
    """
    Test that multiple clients can take turns setting an unlocked resource.
    """
    server, client1 = server_client_pair
    with MultiplayerClient(port=server.port) as client2:
        assert client1.try_set_resource_value("scene", scene)
        assert client2.try_set_resource_value("scene", scene)
        assert client1.try_set_resource_value("scene", scene)
        assert client2.try_set_resource_value("scene", scene)

def test_can_set_own_locked(server_client_pair, scene):
    """
    Test that you can set a resource value if you have locked it.
    """
    server, client = server_client_pair
    client.try_lock_resource("scene")
    assert client.try_set_resource_value("scene", scene)


def test_set_value_updates_server_values(server_client_pair, scene):
    """
    Test that setting a resource value updates the server's internal resource
    map.
    """
    server, client = server_client_pair
    client.try_set_resource_value("scene", scene)
    server_scene = server.multiplayer_service.resources.get("scene")
    assert str(scene) == str(server_scene)


def test_subscribe_value_update(server_client_pair):
    """
    Test that resource value updates can be subscribed.
    """
    server, client = server_client_pair
    client.subscribe_all_value_updates()


def test_set_value_sends_update(server_client_pair, scene):
    """
    Test that setting a resource value is propagated back to the client.
    """
    server, client = server_client_pair
    client.subscribe_all_value_updates()
    time.sleep(0.1)
    client.try_set_resource_value("scene", scene)
    time.sleep(0.1)
    recv_scene = client.resources.get("scene")
    assert str(scene) == str(recv_scene)


def test_server_sends_initial_values(server_client_pair, scene):
    """
    Test that subscribing resource values sends any resources values that have
    already been set.
    """
    server, client = server_client_pair
    client.try_set_resource_value("scene", scene)
    time.sleep(0.1)
    assert client.resources.get("scene") is None
    client.subscribe_all_value_updates()
    time.sleep(0.1)
    assert str(client.resources.get("scene")) == str(scene)


def test_cant_lock_other_locked(server_client_pair):
    """
    Test that you cannot lock a resource that is locked by someone else.
    """
    server, client1 = server_client_pair
    with MultiplayerClient(port=server.port) as client2:
        client1.join_multiplayer("main")
        client2.join_multiplayer("other")
        client2.try_lock_resource("scene")
        assert not client1.try_lock_resource("scene")


def test_cant_release_other_lock(server_client_pair):
    """
    Test that you cannot release a resource that is locked by someone else.
    """
    server, client1 = server_client_pair
    with MultiplayerClient(port=server.port) as client2:
        client1.join_multiplayer("main")
        client2.join_multiplayer("other")
        client2.try_lock_resource("scene")
        assert not client1.try_release_resource("scene")


def test_cant_set_other_locked(server_client_pair, scene):
    """
    Test that you cannot set a resource that is locked by someone else.
    """
    server, client1 = server_client_pair
    with MultiplayerClient(port=server.port) as client2:
        client1.join_multiplayer("main")
        client2.join_multiplayer("other")
        client2.try_lock_resource("scene")
        assert not client1.try_set_resource_value("scene", scene)


def test_cant_set_non_value(server_client_pair):
    """
    Test that setting a value to a non-grpc value raises the appropriate
    :exception.
    """
    server, client = server_client_pair
    with pytest.raises(TypeError):
        client.try_set_resource_value("scene", "hello")
