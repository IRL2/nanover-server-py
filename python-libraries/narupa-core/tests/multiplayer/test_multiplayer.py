# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Integration tests of the multiplayer server with the reference multiplayer client.
"""
import threading
import time

import pytest

from narupa.multiplayer.multiplayer_client import MultiplayerClient
from narupa.multiplayer.multiplayer_server import MultiplayerServer

CONNECT_WAIT_TIME = 0.01
IMMEDIATE_REPLY_WAIT_TIME = 0.01
CLIENT_SEND_INTERVAL = 1 / 60


@pytest.fixture
def server_client_pair():
    """
    Provides a multiplayer server hosting on an available port on localhost,
    and a multiplayer client connected to it.
    """
    server = MultiplayerServer(address='localhost', port=0)
    client = MultiplayerClient.insecure_channel(
        port=server.port,
        send_interval=CLIENT_SEND_INTERVAL
    )

    with client, server:
        yield server, client


@pytest.fixture
def scene():
    """
    Provides scene test data.
    """
    return {
        'position': {"x": 1., "y": 1., "z": 1.},
        'rotation': {"x": 0., "y": 0., "z": 0., "w": 1.},
        'scale': 1.,
    }


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
    with MultiplayerClient.insecure_channel(port=server.port) as client2:
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
    server_scene = server._multiplayer_service.copy_state().get("scene")
    assert scene == server_scene


def test_remove_key_updates_server_keys(server_client_pair, scene):
    """
    Test that setting a resource value updates the server's internal resource
    map.
    """
    server, client = server_client_pair
    client.try_set_resource_value("scene", scene)
    server_scene = server._multiplayer_service.copy_state().get("scene")
    assert scene == server_scene
    client.try_remove_resource_key("scene")
    server_scene = server._multiplayer_service.copy_state().get("scene")
    assert server_scene is None


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
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    client.try_set_resource_value("scene", scene)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    recv_scene = client.resources.get("scene")
    assert scene == recv_scene


def test_remove_key_sends_update(server_client_pair, scene):
    """
    Test that removing a resource key is propagated back to the client.
    """
    server, client = server_client_pair
    client.subscribe_all_value_updates()
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    client.try_set_resource_value("scene", scene)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    recv_scene = client.resources.get("scene")
    assert scene == recv_scene
    client.try_remove_resource_key("scene")
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    recv_scene = client.resources.get("scene")
    assert recv_scene is None


def test_server_sends_initial_values(server_client_pair, scene):
    """
    Test that subscribing resource values sends any resources values that have
    already been set.
    """
    server, client = server_client_pair
    client.try_set_resource_value("scene", scene)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    assert client.resources.get("scene") is None
    client.subscribe_all_value_updates()
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    assert client.resources.get("scene") == scene


def test_cant_lock_other_locked(server_client_pair):
    """
    Test that you cannot lock a resource that is locked by someone else.
    """
    server, client1 = server_client_pair
    with MultiplayerClient.insecure_channel(port=server.port) as client2:
        client1.join_multiplayer("main")
        client2.join_multiplayer("other")
        client1._player_id = '1'
        client2._player_id = '2'
        client2.try_lock_resource("scene")
        assert not client1.try_lock_resource("scene")


def test_cant_release_other_lock(server_client_pair):
    """
    Test that attempting to release a lock you no longer hold but is now held
    by someone else will succeed but leave the locks untouched.
    """
    server, client1 = server_client_pair
    with MultiplayerClient.insecure_channel(port=server.port) as client2:
        client1.join_multiplayer("main")
        client2.join_multiplayer("other")
        assert client2.try_lock_resource("scene")
        assert client1.try_release_resource("scene")
        assert not client1.try_lock_resource("scene")


def test_cant_set_other_locked(server_client_pair, scene):
    """
    Test that you cannot set a resource that is locked by someone else.
    """
    server, client1 = server_client_pair
    with MultiplayerClient.insecure_channel(port=server.port) as client2:
        client1.join_multiplayer("main")
        client2.join_multiplayer("other")
        client2.try_lock_resource("scene")
        assert not client1.try_set_resource_value("scene", scene)


def test_cant_set_non_basic_type(server_client_pair):
    """
    Test that setting a value to a non-basic type raises a ValueError.
    """
    server, client = server_client_pair

    class UserType:
        pass

    with pytest.raises(ValueError):
        client.try_set_resource_value("scene", UserType())


@pytest.mark.parametrize('update_interval', (1 / 10, 1 / 30, 1 / 60))
def test_subscribe_value_sends_initial_immediately(server_client_pair,
                                                   update_interval):
    """
    Test that subscribing values before any have been sent will immediately
    send the first update regardless of interval.
    """
    server, client = server_client_pair
    client.join_multiplayer("main", join_streams=False)
    client.subscribe_all_value_updates(interval=update_interval)

    test_value = "hello"

    time.sleep(CONNECT_WAIT_TIME)
    client.try_set_resource_value("test", test_value)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    assert client.resources.get("test") == test_value


@pytest.mark.parametrize('update_interval', (.5, .2, .1))
def test_subscribe_value_interval(server_client_pair, update_interval):
    """
    Test that value updates are sent at the requested interval.
    """
    server, client = server_client_pair
    client.join_multiplayer("main", join_streams=False)
    client.subscribe_all_value_updates(interval=update_interval)
    test_values = [f"hello {i}" for i in range(5)]
    time.sleep(CONNECT_WAIT_TIME)

    client.try_set_resource_value("test", test_values[0])
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    client.try_set_resource_value("test", test_values[1])
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    assert client.resources.get("test") == test_values[0]

    client.try_set_resource_value("test", test_values[2])
    time.sleep(update_interval)
    assert client.resources.get("test") == test_values[2]


@pytest.mark.parametrize('lock_duration', (.5, 1, 2))
def test_lock_durations(server_client_pair, lock_duration):
    """
    Test that locks expire roughly after the requested duration has passed.
    """
    server, client1 = server_client_pair
    with MultiplayerClient.insecure_channel(port=server.port) as client2:
        client1.join_multiplayer("main")
        client2.join_multiplayer("other")

        client1.try_lock_resource("test", duration=lock_duration)
        time.sleep(lock_duration * .9)
        assert not client2.try_lock_resource("test")
        time.sleep(lock_duration * .2)
        assert client2.try_lock_resource("test")


@pytest.mark.parametrize('lock_duration', (.5, 1, 2))
def test_lock_duration_extend(server_client_pair, lock_duration):
    """
    Test that relocking a key updates the lock duration.
    """
    server, client1 = server_client_pair
    with MultiplayerClient.insecure_channel(port=server.port) as client2:
        client1.join_multiplayer("main")
        client2.join_multiplayer("other")

        client1.try_lock_resource("test", duration=lock_duration)

        # relock halfway through the existing lock
        time.sleep(lock_duration * .5)
        client1.try_lock_resource("test", duration=lock_duration)
        # check lock remains after the initially requested duration expires
        time.sleep(lock_duration * .6)
        assert not client2.try_lock_resource("test")
        # check lock then expires after second requested duration expires
        time.sleep(lock_duration * .5)
        assert client2.try_lock_resource("test")


@pytest.mark.timeout(3)
def test_repeated_disconnect_frees_workers():
    """
    Test that disconnecting frees workers on the server. With just enough
    workers to service one subscription and one request, the server will get
    stuck queueing requests if any subscriptions are left hanging.
    """
    # 2 workers: 1 for value subscription, 1 for single request
    with MultiplayerServer(address='localhost', port=0, max_workers=2) as server:
        for _ in range(32):
            with MultiplayerClient.insecure_channel(port=server.port) as client:
                client.join_multiplayer("test", join_streams=False)
                time.sleep(CONNECT_WAIT_TIME)
                client.subscribe_all_value_updates()
                client.try_set_resource_value("test", 0)


@pytest.mark.timeout(3)
def test_repeated_connections_stall_server():
    """
    Test that too many clients not disconnecting will cause a queue in the
    server's executor.
    """
    with MultiplayerServer(address='localhost', port=0, max_workers=2) as server:
        def loads_of_clients():
            clients = []
            try:
                for _ in range(32):
                    client = MultiplayerClient.insecure_channel(port=server.port)
                    clients.append(client)
                    client.join_multiplayer("test")
            finally:
                for client in clients:
                    client.close()

        threading.Thread(target=loads_of_clients, daemon=True).start()
        time.sleep(CONNECT_WAIT_TIME)
        assert server.server._state.thread_pool._work_queue.qsize() > 0
