import time

import grpc
import pytest
from google.protobuf.struct_pb2 import Value

from narupa.multiplayer import MultiplayerServer
from .test_frame_server import simple_frame_data, frame_server
from .imd.test_imd_server import imd_server, interaction
from narupa.app.client import NarupaImdClient
import numpy as np

CLIENT_WAIT_TIME = 0.2

TEST_KEY = 'test'
TEST_VALUE = Value(string_value='hi')


@pytest.fixture
def multiplayer_server():
    server = MultiplayerServer(address='localhost', port=0)
    yield server
    server.close()


@pytest.fixture
def client_server(frame_server, imd_server, multiplayer_server):
    with NarupaImdClient(trajectory_port=frame_server.port,
                         imd_port=imd_server.port,
                         multiplayer_port=multiplayer_server.port) as client:
        yield client, frame_server, imd_server, multiplayer_server


def test_receive_frames(client_server, simple_frame_data):
    client, frame_server, imd_server, multiplayer_server = client_server
    time.sleep(CLIENT_WAIT_TIME)
    frame_server.send_frame(0, simple_frame_data)
    time.sleep(CLIENT_WAIT_TIME)
    assert client.latest_frame is not None
    assert client.first_frame is not None
    assert client.latest_frame == client.first_frame
    assert len(client.frames) == 1


def test_receive_multiple_frames(client_server, simple_frame_data):
    client, frame_server, imd_server, multiplayer_server = client_server
    # wait for client to connect.
    time.sleep(CLIENT_WAIT_TIME)
    frame_server.send_frame(0, simple_frame_data)
    frame_server.send_frame(1, simple_frame_data)
    time.sleep(CLIENT_WAIT_TIME)
    assert len(client.frames) == 2


def test_latest_frame_only(frame_server, simple_frame_data):
    """
    Tests that running with the all_frames flag set to false, the client only receives the latest
    frame.
    """
    with NarupaImdClient(all_frames=False, trajectory_port=frame_server.port) as client:
        num_frames = 2
        for i in range(num_frames):
            frame_server.send_frame(0, simple_frame_data)
        time.sleep(CLIENT_WAIT_TIME)
        assert client.latest_frame is not None
        assert client.first_frame is not None
        assert client.latest_frame == client.first_frame
        assert len(client.frames) < num_frames


def test_reconnect_receive(client_server, simple_frame_data):
    client, frame_server, imd_server, multiplayer_server = client_server
    frame_server.send_frame(0, simple_frame_data)
    client.close()
    assert client.latest_frame is None
    client.connect(trajectory_port=frame_server.port, imd_port=imd_server.port)
    time.sleep(CLIENT_WAIT_TIME)
    assert client.latest_frame is not None


def test_close_interaction(client_server, interaction):
    client, frame_server, imd_server, multiplayer_server = client_server
    id = client.start_interaction(interaction)
    time.sleep(CLIENT_WAIT_TIME)
    assert len(imd_server.service.active_interactions) == 1
    client.close()
    time.sleep(CLIENT_WAIT_TIME)
    assert len(imd_server.service.active_interactions) == 0


def test_stop_interaction(client_server, interaction):
    """
    tests that stopping an interaction produces the expected results on the server.
    """
    client, frame_server, imd_server, multiplayer_server = client_server
    id = client.start_interaction(interaction)
    time.sleep(CLIENT_WAIT_TIME)
    assert len(imd_server.service.active_interactions) == 1
    client.stop_interaction(id)
    time.sleep(CLIENT_WAIT_TIME)
    assert len(imd_server.service.active_interactions) == 0


def test_start_interaction(client_server, interaction):
    """
    tests that starting an interaction produces expected results on the server.
    """
    client, frame_server, imd_server, multiplayer_server = client_server
    client.start_interaction(interaction)
    time.sleep(CLIENT_WAIT_TIME)
    assert len(imd_server.service.active_interactions) == 1


def test_update_interaction(client_server, interaction):
    """
    tests updating an interaction produces expected results on the server.
    """
    client, frame_server, imd_server, multiplayer_server = client_server
    id = client.start_interaction(interaction)
    interaction.position = [2, 2, 2]
    client.update_interaction(id, interaction)
    time.sleep(CLIENT_WAIT_TIME)
    assert len(imd_server.service.active_interactions) == 1
    assert np.allclose(list(imd_server.service.active_interactions.values())[0].position, (2, 2, 2))


def test_no_imd(frame_server, multiplayer_server, interaction):
    """
    tests that running an interaction raises an exception.

    An error will only be thrown when the interaction is stopped.
    """
    client = NarupaImdClient(trajectory_port=frame_server.port, multiplayer_port=multiplayer_server.port)
    id = client.start_interaction(interaction)
    with pytest.raises(grpc.RpcError):
        client.stop_interaction(id)
    client.close()


def test_set_multiplayer_value(client_server):
    """
    tests that setting multiplayer value works correctly.
    """
    client, frame_server, imd_server, multiplayer_server = client_server
    client.join_multiplayer("1")

    client.set_shared_value(TEST_KEY, TEST_VALUE)
    time.sleep(CLIENT_WAIT_TIME)
    assert client.latest_multiplayer_values == {TEST_KEY: TEST_VALUE}


def test_set_multiplayer_value_disconnected(frame_server):
    """
    tests that setting multiplayer value when not connected raises expected exception.
    """
    with NarupaImdClient(trajectory_port=frame_server.port) as client:
        with pytest.raises(grpc.RpcError):
            client.set_shared_value(TEST_KEY, TEST_VALUE)


def test_join_multiplayer_disconnected(frame_server):
    """
    Tests that joining multiplayer when not connected raises expected exception.
    """
    with NarupaImdClient(trajectory_port=frame_server.port) as client:
        with pytest.raises(grpc.RpcError):
            client.join_multiplayer("1")


def test_get_shared_resources_disconnected(frame_server):
    """
    Tests that getting shared resources produces empty dictionary.
    """
    # TODO handle lack of connection by throwing an exception.
    with NarupaImdClient(trajectory_port=frame_server.port) as client:
        assert client.latest_multiplayer_values == {}
