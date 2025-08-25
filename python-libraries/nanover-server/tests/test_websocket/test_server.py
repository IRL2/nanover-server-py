from contextlib import contextmanager

import msgpack
import pytest
from websockets.sync.server import Server
from websockets.sync.client import connect

from nanover.app import NanoverImdApplication
from nanover.trajectory import FrameData
from nanover.trajectory.frame_data import PARTICLE_COUNT
from nanover.utilities.change_buffers import DictionaryChange
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
    with connect(f"ws://localhost:{port}") as websocket:

        def send_message(message: dict):
            websocket.send(msgpack.packb(message))

        def recv_message():
            return msgpack.unpackb(websocket.recv())

        yield send_message, recv_message


TEST_FRAME = FrameData()
TEST_FRAME.particle_count = 42

TEST_CHANGE = DictionaryChange(
    updates={"baby.yoda": 2000},
    removals={},
)


@pytest.mark.parametrize("frame", (TEST_FRAME,))
def test_websocket_sends_frame(frame):
    with make_websocket_server() as (app_server, ws):
        with connect_client_to_server(ws) as (send, recv):
            _ = recv()

            app_server._frame_publisher.send_frame(frame_index=1, frame=frame)

            message = recv()
            assert "frame" in message

            frame = message["frame"]
            assert frame[PARTICLE_COUNT] == TEST_FRAME.values[PARTICLE_COUNT]


@pytest.mark.parametrize("change", (TEST_CHANGE,))
def test_websocket_sends_state(change):
    with make_websocket_server() as (app_server, ws):
        with connect_client_to_server(ws) as (send, recv):
            _ = recv()

            app_server.server._state_service.state_dictionary.update_state(None, change)

            message = recv()
            assert "state" in message

            update = message["state"]
            assert update["updates"]["baby.yoda"] == TEST_CHANGE.updates["baby.yoda"]
