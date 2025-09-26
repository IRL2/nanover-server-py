from contextlib import contextmanager
from dataclasses import dataclass

from nanover.app import NanoverImdApplication
from nanover.app.types import AppServer
from nanover.testing import assert_equal_soon
from nanover.testing.utilities import simplify_numpy
from nanover.trajectory import FrameData2
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket import NanoverImdClient
from nanover.websocket.server import WebSocketServer


@contextmanager
def make_connected_server_client_setup():
    with make_websocket_server() as (app_server, ws):
        with connect_client_to_server(ws) as client:
            yield ServerClientSetup(
                app=app_server,
                server=ws,
                client=client,
            )


@contextmanager
def make_app_server():
    with NanoverImdApplication.basic_server(port=0) as app_server:
        yield app_server


@contextmanager
def make_websocket_server():
    with make_app_server() as app_server:
        yield app_server, app_server._server_ws


@contextmanager
def connect_client_to_server(websocket_server: WebSocketServer):
    with NanoverImdClient.from_url(
        f"ws://localhost:{websocket_server.ws_port}"
    ) as client:
        yield client


@dataclass(kw_only=True)
class ServerClientSetup:
    app: AppServer
    server: WebSocketServer
    client: NanoverImdClient

    def assert_frames_synced_soon(self, **kwargs):
        __tracebackhide__ = True  # hide this function in the test traceback
        assert_equal_soon(
            lambda: simplify_numpy(self.client_current_frame.frame_dict),
            lambda: simplify_numpy(self.server_current_frame.frame_dict),
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
        self.server_publish_frame(frame=FrameData2(), frame_index=0)

    def server_publish_frame(self, frame: FrameData2, frame_index: int):
        self.app.frame_publisher.send_frame(frame=frame, frame_index=frame_index)

    @property
    def client_current_frame(self):
        return self.client.current_frame

    @property
    def server_current_frame(self):
        return self.app.frame_publisher.last_frame

    def server_update_state(self, change: DictionaryChange):
        self.app.state_dictionary.update_state(None, change)

    @property
    def server_current_state(self):
        return self.app.state_dictionary.copy_content()

    @property
    def client_current_state(self):
        return self.client._state_dictionary.copy_content()
