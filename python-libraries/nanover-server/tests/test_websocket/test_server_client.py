import pytest
from hypothesis import given, example, strategies as st
from mock import Mock
import ssl

from nanover.websocket.server import WebSocketServer
from nanover.app.imd_app import NanoverImdApplication
from nanover.testing.servers import (
    make_connected_server_client_setup,
    connect_client_to_server,
)
from nanover.testing.utilities import simplify_numpy
from nanover.trajectory import FrameData
from nanover.utilities.change_buffers import DictionaryChange

from nanover.testing import assert_equal_soon
from nanover.testing.strategies import (
    command_arguments,
    state_updates,
    frames,
)


@example(port=0, ssl_=True)
@example(port=80, ssl_=True).xfail(raises=ValueError)
@example(port=70000, ssl_=True).xfail(raises=ValueError)
@given(st.integers(1024, 65535), st.booleans())
def test_websocket_server_instantiation(port, ssl_):
    ssl_context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER) if ssl_ else None

    with WebSocketServer.basic_server(
        NanoverImdApplication(), port=port, ssl=ssl_context
    ) as ws_server:
        assert ws_server.ws_port is not None
        port = ws_server.ws_port
        if ssl_:
            assert ws_server.wss_port is not None

        with WebSocketServer.basic_server(
            NanoverImdApplication(), port=port, ssl=ssl_context
        ) as second_server:
            assert second_server.ws_port is not None
            assert second_server.ws_port != ws_server.ws_port

            if ssl_:
                assert second_server.wss_port is not None


@pytest.fixture(scope="module")
def reusable_setup():
    def echo(**arguments):
        return arguments

    with make_connected_server_client_setup() as setup:
        setup.app.register_command("test/identity", echo)
        yield setup


@pytest.fixture(scope="module")
def reusable_setup_two_clients():
    with make_connected_server_client_setup() as setup:
        with connect_client_to_server(setup.server) as client:
            yield setup, client


@given(frame=frames())
def test_websocket_sends_frame(reusable_setup, frame):
    reusable_setup.server_publish_frame(frame)
    reusable_setup.assert_frames_synced_soon()


@given(frame=frames())
def test_websocket_sends_frame_two_clients(reusable_setup_two_clients, frame):
    reusable_setup, client2 = reusable_setup_two_clients

    reusable_setup.server_publish_frame(frame)
    reusable_setup.assert_frames_synced_soon()

    assert_equal_soon(
        lambda: simplify_numpy(reusable_setup.client.current_frame.frame_dict),
        lambda: simplify_numpy(client2.current_frame.frame_dict),
    )

    reusable_setup.server_publish_frame_reset()
    reusable_setup.assert_frames_synced_soon()

    assert_equal_soon(
        lambda: simplify_numpy(reusable_setup.client.current_frame.frame_dict),
        lambda: simplify_numpy(client2.current_frame.frame_dict),
    )


@given(updates=state_updates())
def test_websocket_sends_state(reusable_setup, updates):
    """
    Test that state updates made directly on the server are accurately reflected on the client.
    """
    change = DictionaryChange(updates=updates)
    reusable_setup.server_update_state(change)

    assert_equal_soon(
        lambda: pick(reusable_setup.server_current_state, updates),
        lambda: updates,
    )

    reusable_setup.assert_states_synced_soon()


@given(arguments=command_arguments())
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


@given(frame=frames())
def test_client_frame_reset(reusable_setup, frame):
    reusable_setup.server_publish_frame(frame)
    reusable_setup.assert_frames_synced_soon()

    reusable_setup.server_publish_frame_reset()
    reusable_setup.assert_frames_synced_soon()


def pick(dictionary, keys):
    return {key: value for key, value in dictionary.items() if key in keys}
