import time
import pytest
from mock import Mock

from nanover.app import NanoverImdApplication
from nanover.essd.utils import get_broadcastable_test_ip
from nanover.imd import ParticleInteraction
from nanover.testing import assert_equal_soon, assert_in_soon
from nanover.trajectory import FrameData, keys

from .test_frame_server import simple_frame_data, disjoint_frame_data
from nanover.websocket import NanoverImdClient
import numpy as np

TEST_KEY = "test"
TEST_VALUE = "hi"


@pytest.fixture
def interaction():
    return ParticleInteraction()


@pytest.fixture
def default_args():
    return {"a": 2, "b": [1, 3, 4], "c": True}


@pytest.fixture
def mock_callback(default_args):
    return Mock(return_value=default_args)


@pytest.fixture
def client_server():
    with NanoverImdApplication.basic_server(
        address=get_broadcastable_test_ip(), port=0
    ) as app_server:
        with NanoverImdClient.from_app_server(app_server) as client:
            yield client, app_server


def test_receive_multiple_frames(client_server, simple_frame_data):
    client, app_server = client_server
    app_server.frame_publisher.send_frame(simple_frame_data)
    time.sleep(0.1)
    app_server.frame_publisher.send_frame(simple_frame_data)
    assert_equal_soon(
        lambda: len(client.frames),
        lambda: 2,
    )


def test_current_frame_does_merge(client_server):
    client, app_server = client_server

    first_frame = FrameData()
    first_frame["indices"] = [0, 1, 3]
    first_frame["string"] = "str"

    second_frame = FrameData()
    second_frame["indices"] = [4, 6, 8]
    second_frame["bool"] = False

    app_server.frame_publisher.send_frame(first_frame)
    app_server.frame_publisher.send_frame(second_frame)

    assert_equal_soon(
        lambda: (
            client.current_frame["indices"],
            client.current_frame["bool"],
            client.current_frame["string"],
        ),
        lambda: (
            second_frame["indices"],
            second_frame["bool"],
            first_frame["string"],
        ),
    )


def test_frame_reset(client_server, simple_frame_data, disjoint_frame_data):
    client, app_server = client_server

    app_server.frame_publisher.send_frame(simple_frame_data)

    assert_in_soon(
        lambda: "bool",
        lambda: client.current_frame,
    )

    app_server.frame_publisher.send_clear()
    app_server.frame_publisher.send_frame(disjoint_frame_data)

    assert_in_soon(
        lambda: "number",
        lambda: client.current_frame,
    )
    assert "bool" not in client.current_frame


def test_close_interaction(client_server, interaction):
    client, app_server = client_server
    client.start_interaction(interaction)

    assert_equal_soon(
        lambda: len(app_server.imd.active_interactions),
        lambda: 1,
    )

    client.close()

    assert_equal_soon(
        lambda: len(app_server.imd.active_interactions),
        lambda: 0,
    )


def test_stop_interaction(client_server, interaction):
    """
    tests that stopping an interaction produces the expected results on the server.
    """
    client, app_server = client_server
    id = client.start_interaction(interaction)

    assert_equal_soon(
        lambda: len(app_server.imd.active_interactions),
        lambda: 1,
    )

    client.stop_interaction(id)

    assert_equal_soon(
        lambda: len(app_server.imd.active_interactions),
        lambda: 0,
    )


def test_start_interaction(client_server, interaction):
    """
    tests that starting an interaction produces expected results on the server.
    """
    client, app_server = client_server
    client.start_interaction(interaction)

    assert_equal_soon(
        lambda: len(app_server.imd.active_interactions),
        lambda: 1,
    )


def test_update_interaction(client_server, interaction):
    """
    tests updating an interaction produces expected results on the server.
    """
    client, app_server = client_server
    id = client.start_interaction(interaction)
    interaction.scale = 0
    interaction.position = [2, 2, 2]
    client.update_interaction(id, interaction)

    assert_equal_soon(
        lambda: list(app_server.imd.active_interactions.values())[0].scale,
        lambda: 0,
    )

    assert np.allclose(
        list(app_server.imd.active_interactions.values())[0].position,
        (2, 2, 2),
    )


def test_set_multiplayer_value(client_server):
    """
    tests that setting multiplayer value works correctly.
    """
    client, app_server = client_server

    client.set_shared_value(TEST_KEY, TEST_VALUE)

    assert_equal_soon(
        lambda: client.latest_multiplayer_values[TEST_KEY],
        lambda: TEST_VALUE,
    )


def test_run_play(client_server, mock_callback):
    client, app_server = client_server
    app_server.register_command(keys.PLAY_COMMAND, mock_callback)
    client.run_play()
    mock_callback.assert_called_once()


def test_run_reset(client_server, mock_callback):
    client, app_server = client_server
    app_server.register_command(keys.RESET_COMMAND, mock_callback)
    client.run_reset()
    mock_callback.assert_called_once()


def test_run_step(client_server, mock_callback):
    client, app_server = client_server
    app_server.register_command(keys.STEP_COMMAND, mock_callback)
    client.run_step()
    mock_callback.assert_called_once()


def test_run_pause(client_server, mock_callback):
    client, app_server = client_server
    app_server.register_command(keys.PAUSE_COMMAND, mock_callback)
    client.run_pause()
    mock_callback.assert_called_once()


frame_str = "frame"
imd_str = "imd"
multiplayer_str = "multiplayer"


def test_available_commands(client_server, mock_callback):
    client, app_server = client_server
    app_server.register_command(frame_str, mock_callback)
    app_server.register_command(imd_str, mock_callback)
    app_server.register_command(multiplayer_str, mock_callback)

    commands = client.update_available_commands()

    assert set(commands.keys()).issuperset({frame_str, imd_str, multiplayer_str})


def test_available_commands_frame_server_only(client_server, mock_callback):
    """
    tests that if other servers are not set up, requesting the available commands
    still works.
    """
    client, app_server = client_server
    app_server.register_command(frame_str, mock_callback)
    commands = client.update_available_commands()
    assert set(commands.keys()).issuperset({frame_str})


def run_client_server_command_test(client, server):
    return_mock = {frame_str: frame_str}
    mock = Mock(return_value=return_mock)
    server.register_command(frame_str, mock)
    result = client.run_command_blocking(frame_str)
    assert result == {frame_str: frame_str}


def test_unknown_command(client_server):
    client, app_server = client_server
    with pytest.raises(RuntimeError):
        client.run_command_blocking("unknown")
