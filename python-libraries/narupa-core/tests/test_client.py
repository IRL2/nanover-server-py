import time

import grpc
import pytest
from mock import Mock

from narupa.multiplayer import MultiplayerServer
from narupa.trajectory.frame_server import PLAY_COMMAND_KEY, RESET_COMMAND_KEY, STEP_COMMAND_KEY, PAUSE_COMMAND_KEY

from .test_frame_server import simple_frame_data, frame_server
from .imd.test_imd_server import imd_server, interaction
from .core.test_grpc_client_server_commands import mock_callback, default_args
from narupa.app.client import NarupaImdClient
import numpy as np

CLIENT_WAIT_TIME = 0.2

TEST_KEY = 'test'
TEST_VALUE = 'hi'


@pytest.fixture
def multiplayer_server():
    with MultiplayerServer(address='localhost', port=0) as server:
        yield server


@pytest.fixture
def client_server(frame_server, imd_server, multiplayer_server):
    with NarupaImdClient(trajectory_port=frame_server.port,
                         imd_port=imd_server.port,
                         multiplayer_port=multiplayer_server.port) as client:
        yield client, frame_server, imd_server, multiplayer_server


@pytest.fixture
def client_frame_server(frame_server):
    with NarupaImdClient(trajectory_port=frame_server.port) as client:
        yield client, frame_server


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


def test_run_play(client_frame_server, mock_callback):
    client, frame_server = client_frame_server
    frame_server.register_command(PLAY_COMMAND_KEY, mock_callback)
    client.run_play()
    mock_callback.assert_called_once()


def test_run_reset(client_frame_server, mock_callback):
    client, frame_server = client_frame_server
    frame_server.register_command(RESET_COMMAND_KEY, mock_callback)
    client.run_reset()
    mock_callback.assert_called_once()


def test_run_step(client_frame_server, mock_callback):
    client, frame_server = client_frame_server
    frame_server.register_command(STEP_COMMAND_KEY, mock_callback)
    client.run_step()
    mock_callback.assert_called_once()


def test_run_pause(client_frame_server, mock_callback):
    client, frame_server = client_frame_server
    frame_server.register_command(PAUSE_COMMAND_KEY, mock_callback)
    client.run_pause()
    mock_callback.assert_called_once()


frame_str = "frame"
imd_str = "imd"
multiplayer_str = "multiplayer"


def test_available_commands(client_server, mock_callback):
    client, frame_server, imd_server, multiplayer_server = client_server
    frame_server.register_command(frame_str, mock_callback)
    imd_server.register_command(imd_str, mock_callback)
    multiplayer_server.register_command(multiplayer_str, mock_callback)

    commands = client.update_available_commands()

    assert len(commands) == 3
    assert set(commands.keys()) == {frame_str, imd_str, multiplayer_str}


def test_available_commands_frame_server_only(client_frame_server, mock_callback):
    """
    tests that if other servers are not set up, requesting the available commands
    still works.
    """
    client, frame_server = client_frame_server
    frame_server.register_command(frame_str, mock_callback)
    commands = client.update_available_commands()
    assert len(commands) == 1
    assert set(commands.keys()) == {frame_str}


def run_client_server_command_test(client, server):
    return_mock = {frame_str: frame_str}
    mock = Mock(return_value=return_mock)
    server.register_command(frame_str, mock)
    client.update_available_commands()
    result = client.run_command(frame_str)
    assert result == {frame_str: frame_str}


def test_run_frame_command_generic(client_server):
    """
    tests that the client can run command on the frame server,
    without having to know which server the command needs to go to.
    """
    client, frame_server, imd_server, multiplayer_server = client_server
    run_client_server_command_test(client, frame_server)


def test_run_multiplayer_command_generic(client_server):
    """
    tests that the client can run command on the multiplayer server,
    without having to know which server the command needs to go to.
    """
    client, frame_server, imd_server, multiplayer_server = client_server
    run_client_server_command_test(client, multiplayer_server)


def test_run_imd_command_generic(client_server):
    """
    tests that the client can run command on the imd server,
    without having to know which server the command needs to go to.
    """
    client, frame_server, imd_server, multiplayer_server = client_server
    run_client_server_command_test(client, imd_server)


def test_run_command_multiple_servers(client_server):
    """
    Tests that commands can be generically run on both the frame server and imd server.
    """
    client, frame_server, imd_server, multiplayer_server = client_server

    mock_frame = Mock(return_value={frame_str: frame_str})
    mock_imd = Mock(return_value={imd_str: imd_str})
    frame_server.register_command(frame_str, mock_frame)
    imd_server.register_command(imd_str, mock_imd)

    commands = client.update_available_commands()

    for command in commands:
        client.run_command(command)

    mock_frame.assert_called_once()
    mock_imd.assert_called_once()


def test_unknown_command(client_server):
    client, frame_server, imd_server, multiplayer_server = client_server
    client.update_available_commands()
    with pytest.raises(KeyError):
        client.run_command("unknown")
