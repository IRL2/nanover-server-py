from contextlib import contextmanager
from dataclasses import dataclass

import pytest
from hypothesis import given, strategies as st
from mock import Mock

from nanover.app import NanoverImdApplication
from nanover.trajectory import FrameData
from nanover.trajectory.frame_data import PARTICLE_COUNT
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.client import WebsocketClient
from nanover.websocket.convert import pack_grpc_frame, unpack_dict_frame
from nanover.websocket.server import WebSocketServer

from nanover.testing import assert_equal_soon
from nanover.testing.strategies import command_arguments, state_updates


@dataclass(kw_only=True)
class ServerClientSetup:
    app: NanoverImdApplication
    server: WebSocketServer
    client: WebsocketClient

    def assert_frames_synced_soon(self, **kwargs):
        __tracebackhide__ = True  # hide this function in the test traceback
        assert_equal_soon(
            lambda: self.client_current_frame,
            lambda: self.server_current_frame,
            **kwargs,
        )

    def assert_states_synced_soon(self, **kwargs):
        __tracebackhide__ = True  # hide this function in the test traceback
        assert_equal_soon(
            lambda: self.client_current_state,
            lambda: self.server_current_state,
            **kwargs,
        )

    def server_publish_frame_reset(self):
        self.server_publish_frame(frame=FrameData(), frame_index=0)

    def server_publish_frame(self, frame: FrameData, frame_index: int):
        self.app._frame_publisher.send_frame(frame=frame, frame_index=frame_index)

    @property
    def client_current_frame(self):
        return self.client.current_frame

    @property
    def server_current_frame(self):
        return unpack_dict_frame(
            pack_grpc_frame(FrameData(self.app._frame_publisher.last_frame))
        )

    def server_update_state(self, change: DictionaryChange):
        self.app.server._state_service.state_dictionary.update_state(None, change)

    @property
    def server_current_state(self):
        return self.app.server._state_service.state_dictionary.copy_content()

    @property
    def client_current_state(self):
        return self.client.state_dictionary.copy_content()


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
    with WebsocketClient.from_url(
        f"ws://localhost:{websocket_server.ws_port}"
    ) as client:
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
    reusable_setup.server_publish_frame(frame_index=frame_index, frame=frame)
    reusable_setup.assert_frames_synced_soon()


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


@given(frame_index=st.integers(min_value=1, max_value=2**32 - 1))
@pytest.mark.parametrize("frame", (TEST_FRAME,))
def test_client_frame_reset(reusable_setup, frame, frame_index):
    reusable_setup.server_publish_frame(frame_index=frame_index, frame=frame)
    reusable_setup.assert_frames_synced_soon()

    reusable_setup.server_publish_frame_reset()
    reusable_setup.assert_frames_synced_soon()


def pick(dictionary, keys):
    return {key: value for key, value in dictionary.items() if key in keys}
