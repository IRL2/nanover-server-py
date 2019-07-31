# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Integration tests of the multiplayer server with the reference multiplayer client.
"""

import time

import grpc
import pytest

from narupa.multiplayer.multiplayer_client import MultiplayerClient
from narupa.multiplayer.multiplayer_server import MultiplayerServer
from narupa.protocol.multiplayer.multiplayer_pb2 import Avatar, AvatarComponent
import narupa.protocol.multiplayer.multiplayer_pb2 as mult_proto


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
def multiplayer_server_channel_send_self(multiplayer_server_send_self):
    channel = grpc.insecure_channel(f'localhost:{multiplayer_server_send_self.port}')

    with channel:
        yield multiplayer_server_send_self, channel


@pytest.fixture
def multiplayer_server_client_send_self(multiplayer_server_send_self):
    client = MultiplayerClient(port=multiplayer_server_send_self.port)
    yield multiplayer_server, client
    client.close()


@pytest.fixture
def avatar():
    components = [AvatarComponent(name="Head", position=[0, 0, 1], rotation=[1, 1, 1, 1])]
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


def test_publish_avatar_not_joined(multiplayer_server_client):
    multiplayer_server, multiplayer_client = multiplayer_server_client
    with pytest.raises(RuntimeError):
        result = multiplayer_client.join_avatar_stream()
        print('Result', result)


def test_join_avatar_stream(multiplayer_server_client):
    multiplayer_server, multiplayer_client = multiplayer_server_client
    multiplayer_client.join_multiplayer(player_name="user", join_streams=False)
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


def test_join_publish_avatar(multiplayer_server_client_send_self):
    multiplayer_server, multiplayer_client = multiplayer_server_client_send_self
    multiplayer_client.join_multiplayer(player_name="user", join_streams=False)
    multiplayer_client.join_avatar_publish()


def test_publish_avatar(multiplayer_server_client_send_self):
    multiplayer_server, multiplayer_client = multiplayer_server_client_send_self
    player_id = multiplayer_client.join_multiplayer(player_name="user", join_streams=True)
    time.sleep(0.05)

    components = [AvatarComponent(name="Head", position=[0, 0, 1], rotation=[1, 1, 1, 1])]
    avatar = Avatar(player_id=player_id, component=components)
    multiplayer_client.publish_avatar(avatar)

    time.sleep(0.05)
    assert len(multiplayer_client.current_avatars) == 1


def test_publish_avatar_multiple_transmission(multiplayer_server_client_send_self, avatar):
    multiplayer_server, multiplayer_client = multiplayer_server_client_send_self
    player_id = multiplayer_client.join_multiplayer(player_name="user", join_streams=True)
    time.sleep(0.05)

    multiplayer_client.publish_avatar(avatar)
    avatar.component[0].position[:] = [0, 0, 2]
    multiplayer_client.publish_avatar(avatar)
    avatar.component[0].position[:] = [0, 0, 3]
    multiplayer_client.publish_avatar(avatar)
    time.sleep(0.05)
    assert len(multiplayer_client.current_avatars) == 1
    assert multiplayer_client.current_avatars[player_id].component[0].position == [0, 0, 3]


@pytest.fixture
def test_scene():
    rotation = [0, 0, 0, 1]
    scale = 1
    properties = None#SceneProperties(position=[1, 1, 1], rotation=rotation, scale=scale)
    return properties


def test_set_scene_properties(multiplayer_server_client_send_self, test_scene):
    multiplayer_server, multiplayer_client = multiplayer_server_client_send_self
    multiplayer_client.join_multiplayer(player_name="user", join_streams=True)
    time.sleep(0.05)

    success = multiplayer_client.set_scene_properties(test_scene)
    assert success


def test_set_scene_properties_received(multiplayer_server_client_send_self, test_scene):
    multiplayer_server, multiplayer_client = multiplayer_server_client_send_self
    multiplayer_client.join_multiplayer(player_name="user", join_streams=True)
    time.sleep(0.05)

    multiplayer_client.set_scene_properties(test_scene)
    time.sleep(0.05)
    received_properties = multiplayer_client.scene_properties
    assert received_properties.properties.position == [1, 1, 1]


def test_self_sending_off_scene_properties(multiplayer_server_client, test_scene):
    """
    Tests that multiplayer client does not receive scene properties it publishes, by default.
    """
    multiplayer_server, multiplayer_client = multiplayer_server_client
    multiplayer_client.join_multiplayer(player_name="user", join_streams=True)
    time.sleep(0.05)

    multiplayer_client.set_scene_properties(test_scene)
    time.sleep(0.05)
    assert multiplayer_client.scene_properties is None


def test_edit_scene_receive_late(multiplayer_server_channel_send_self, test_scene):
    """
    Tests that a second multiplayer client receives the scene property changes if it joins late.

    """
    multiplayer_server, channel = multiplayer_server_channel_send_self
    multiplayer_client = MultiplayerClient(channel=channel)
    multiplayer_client.join_multiplayer(player_name="user", join_streams=True)
    time.sleep(0.05)

    multiplayer_client.set_scene_properties(test_scene)
    time.sleep(0.05)
    assert multiplayer_client.scene_properties is not None

    another_client = MultiplayerClient(channel=channel)
    another_client.join_multiplayer(player_name="someotherguy", join_streams=True)
    time.sleep(0.05)
    assert another_client.scene_properties is not None


def test_edit_scene_locked(multiplayer_server_channel_send_self, test_scene):
    """
    Tests that a second client cannot edit the scene property changes if it is locked

    """
    multiplayer_server, channel = multiplayer_server_channel_send_self
    multiplayer_client = MultiplayerClient(channel=channel)
    multiplayer_client.join_multiplayer(player_name="user", join_streams=True)
    time.sleep(0.05)

    multiplayer_client.try_lock_scene()

    another_client = MultiplayerClient(channel=channel)
    another_client.join_multiplayer(player_name="someotherguy", join_streams=True)
    edit_success = another_client.set_scene_properties(test_scene)
    assert edit_success is False
