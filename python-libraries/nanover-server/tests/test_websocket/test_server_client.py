from contextlib import contextmanager

import msgpack
import pytest
from mock import Mock
from websockets.sync.server import Server
from websockets.sync.client import connect

from nanover.app import NanoverImdApplication
from nanover.testing import assert_equal_soon
from nanover.trajectory import FrameData
from nanover.trajectory.frame_data import PARTICLE_COUNT
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.client import WebsocketClient
from nanover.websocket.server import serve_from_app_server


@contextmanager
def make_app_server():
    with NanoverImdApplication.basic_server(port=0) as app_server:
        yield app_server


@contextmanager
def make_websocket_server():
    with make_app_server() as app_server:
        with serve_from_app_server(app_server) as (wss, ws):
            yield app_server, ws
        ws.shutdown()


@contextmanager
def connect_client_to_server(server: Server):
    port = server.socket.getsockname()[1]
    with WebsocketClient(f"ws://localhost:{port}") as client:
        yield client


@contextmanager
def make_connected_server_client_pair():
    with make_websocket_server() as (app_server, ws):
        with connect_client_to_server(ws) as client:
            yield app_server, client


TEST_FRAME = FrameData()
TEST_FRAME.particle_count = 42

TEST_CHANGE = DictionaryChange(
    updates={"baby.yoda": 2000},
    removals={},
)


TEST_ARGUMENTS = {
    "adult_yoda": 9000,
}


@pytest.mark.parametrize("frame", (TEST_FRAME,))
def test_websocket_sends_frame(frame):
    with make_connected_server_client_pair() as (app_server, client):
        app_server._frame_publisher.send_frame(frame_index=1, frame=frame)
        assert_equal_soon(
            lambda: client.frame.get(PARTICLE_COUNT, None),
            lambda: TEST_FRAME.values[PARTICLE_COUNT],
        )


@pytest.mark.parametrize("change", (TEST_CHANGE,))
def test_websocket_sends_state(change):
    with make_connected_server_client_pair() as (app_server, client):
        service = app_server.server._state_service
        service.state_dictionary.update_state(None, change)
        assert_equal_soon(
            service.state_dictionary.copy_content,
            client.state_dictionary.copy_content,
        )


@pytest.mark.parametrize("arguments", (TEST_ARGUMENTS,))
def test_command(arguments):
    mock_callback = Mock()
    with make_connected_server_client_pair() as (app_server, client):
        app_server.server.register_command("test/command", lambda **kwargs: kwargs)
        client.run_command("test/command", TEST_ARGUMENTS, mock_callback)

        assert_equal_soon(
            lambda: mock_callback.call_args and mock_callback.call_args.args[0],
            lambda: arguments,
        )
