from contextlib import contextmanager
from dataclasses import dataclass

import pytest
from hypothesis import given, strategies as st
from mock import Mock

from nanover.app import NanoverImdApplication
from nanover.testing import assert_equal_soon
from nanover.trajectory import FrameData
from nanover.trajectory.frame_data import PARTICLE_COUNT, FRAME_INDEX
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.client import WebsocketClient
from nanover.websocket.server import WebSocketServer


from data_strategies import packable_structures, dictionary_keys


@dataclass(kw_only=True)
class ServerClientSetup:
    app: NanoverImdApplication
    server: WebSocketServer
    client: WebsocketClient


@pytest.fixture(scope="module")
def reusable_setup():
    with make_connected_server_client_setup() as setup:
        setup.app.server.register_command("test/identity", lambda **kwargs: kwargs)
        yield setup


@contextmanager
def make_app_server():
    with NanoverImdApplication.basic_server(port=0) as app_server:
        yield app_server


@contextmanager
def make_websocket_server():
    with make_app_server() as app_server:
        with WebSocketServer.basic_server(app_server) as websocket_server:
            yield app_server, websocket_server


@contextmanager
def connect_client_to_server(websocket_server: WebSocketServer):
    with WebsocketClient(f"ws://localhost:{websocket_server.ws_port}") as client:
        yield client


@contextmanager
def make_connected_server_client_setup():
    with make_websocket_server() as (app_server, ws):
        with connect_client_to_server(ws) as client:
            yield ServerClientSetup(
                app=app_server,
                server=ws,
                client=client,
            )


TEST_FRAME = FrameData()
TEST_FRAME.particle_count = 42


@given(frame_index=st.integers(min_value=1, max_value=2**32 - 1))
@pytest.mark.parametrize("frame", (TEST_FRAME,))
def test_websocket_sends_frame(reusable_setup, frame, frame_index):
    reusable_setup.app._frame_publisher.send_frame(frame_index=frame_index, frame=frame)
    assert_equal_soon(
        lambda: pick(
            reusable_setup.client.current_frame, {PARTICLE_COUNT, FRAME_INDEX}
        ),
        lambda: {
            PARTICLE_COUNT: TEST_FRAME.values[PARTICLE_COUNT],
            FRAME_INDEX: frame_index,
        },
    )


@pytest.mark.parametrize("frame", (TEST_FRAME,))
def test_websocket_sends_frame_two_clients(frame):
    with make_websocket_server() as (app_server, websocket_server):
        with connect_client_to_server(
            websocket_server
        ) as client1, connect_client_to_server(websocket_server) as client2:
            app_server._frame_publisher.send_frame(frame_index=1, frame=frame)

            assert_equal_soon(
                lambda: (
                    client1.current_frame.get(PARTICLE_COUNT, None),
                    client2.current_frame.get(PARTICLE_COUNT, None),
                ),
                lambda: (TEST_FRAME.values[PARTICLE_COUNT],) * 2,
            )

            app_server._frame_publisher.send_frame(frame_index=0, frame=FrameData())

            assert_equal_soon(
                lambda: client1.current_frame.get(PARTICLE_COUNT, None),
                lambda: client2.current_frame.get(PARTICLE_COUNT, None),
            )


@given(
    updates=st.dictionaries(
        keys=dictionary_keys(), values=packable_structures(), max_size=5
    )
)
def test_websocket_sends_state(reusable_setup, updates):
    """
    Test that state updates made directly on the server are accurately reflected on the client.
    """
    change = DictionaryChange(updates=updates)

    service = reusable_setup.app.server._state_service
    service.state_dictionary.update_state(None, change)

    assert_equal_soon(
        lambda: pick(service.state_dictionary.copy_content(), updates),
        lambda: updates,
    )

    assert_equal_soon(
        lambda: pick(service.state_dictionary.copy_content(), updates),
        lambda: pick(reusable_setup.client.state_dictionary.copy_content(), updates),
    )


@given(
    arguments=st.dictionaries(
        keys=dictionary_keys(), values=packable_structures(), max_size=5
    )
)
def test_echo_command(reusable_setup, arguments):
    """
    Test that preprepared identity command successfully returns unaltered and arbitrary arguments the command
    is called with.
    """
    mock_callback = Mock()
    reusable_setup.client.run_command("test/identity", arguments, mock_callback)

    assert_equal_soon(
        lambda: mock_callback.call_args and mock_callback.call_args.args[0],
        lambda: arguments,
    )


def pick(dictionary, keys):
    return {key: value for key, value in dictionary.items() if key in keys}
