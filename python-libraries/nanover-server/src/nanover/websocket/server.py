from concurrent.futures import ThreadPoolExecutor
from contextlib import suppress, contextmanager
from ssl import SSLContext
from typing import Optional

import msgpack
from websockets import ConnectionClosedOK

from nanover.app import NanoverImdApplication
from nanover.trajectory.frame_data import FrameData
from nanover.utilities.change_buffers import DictionaryChange
from nanover.utilities.cli import CancellationToken
from websockets.sync.server import serve, ServerConnection, Server

from nanover.websocket.convert import convert_frame


@contextmanager
def serve_from_app_server(
    app_server: NanoverImdApplication,
    *,
    ssl: Optional[SSLContext] = None,
    cancellation: Optional[CancellationToken] = None,
):
    cancellation = cancellation or CancellationToken()

    def handle_client(websocket: ServerConnection):
        WebSocketClientHandler(app_server, websocket, cancellation).listen()

    with (
        serve(handle_client, "0.0.0.0", 0, ssl=ssl) as wss,
        serve(handle_client, "0.0.0.0", 0) as ws,
    ):
        cancellation.subscribe_cancellation(wss.shutdown)
        cancellation.subscribe_cancellation(ws.shutdown)

        if app_server.running_discovery:
            hub = app_server._service_hub
            hub.add_service("ws", get_server_port(ws))
            if ssl is not None:
                hub.add_service("wss", get_server_port(wss))
            app_server._update_discovery_services()

        threads = ThreadPoolExecutor(max_workers=2)
        threads.submit(wss.serve_forever)
        threads.submit(ws.serve_forever)
        yield wss, ws


class WebSocketClientHandler:
    def __init__(
        self,
        app_server: NanoverImdApplication,
        websocket: ServerConnection,
        cancellation: CancellationToken,
    ):
        self.app_server = app_server
        self.websocket = websocket

        self.cancellation = CancellationToken()
        self.cancellation.subscribe_cancellation(self.websocket.close)
        cancellation.subscribe_cancellation(self.cancellation.cancel)

    def close(self):
        self.cancellation.cancel()

    @property
    def frame_publisher(self):
        return self.app_server._frame_publisher

    @property
    def state_dictionary(self):
        return self.app_server.server._state_service.state_dictionary

    def run_command(self, name: str, arguments: Optional[dict] = None):
        return self.app_server.server._command_service.run_command(
            name, arguments or {}
        )

    def send_frame(self, frame: FrameData):
        self.send_message({"frame": convert_frame(frame)})

    def send_state_update(self, change: DictionaryChange):
        self.send_message(
            {
                "state": {
                    "updates": change.updates,
                    "removals": list(change.removals),
                },
            }
        )

    def send_message(self, message):
        self.websocket.send(msgpack.packb(message))

    def recv_message(self, message: dict):
        def handle_state_update(update):
            change = DictionaryChange(
                updates=update.get("updates", {}),
                removals=update.get("removals", set()),
            )
            self.state_dictionary.update_state(None, change)

        def handle_command_request(request):
            name, arguments = request.get("name"), request.get("arguments", {})
            result = self.run_command(name, arguments)
            response = {"request": request, "response": result}
            return response

        if "state" in message:
            handle_state_update(message["state"])
        if "command" in message:
            requests = (
                request.get("request", request) for request in message["command"]
            )
            responses = [handle_command_request(request) for request in requests]
            self.send_message({"command": responses})

    def listen(self, frame_interval=1 / 30, state_interval=1 / 30):
        def send_frames():
            for response in self.frame_publisher.subscribe_latest_frames(
                frame_interval=frame_interval,
                cancellation=self.cancellation,
            ):
                self.send_frame(FrameData(response.frame))

        def send_updates():
            with self.state_dictionary.get_change_buffer() as change_buffer:
                self.cancellation.subscribe_cancellation(change_buffer.freeze)
                for change in change_buffer.subscribe_changes(interval=state_interval):
                    self.send_state_update(change)

        def recv_all():
            with suppress(ConnectionClosedOK):
                for data in self.websocket:
                    self.recv_message(msgpack.unpackb(data))
            self.cancellation.cancel()

        threads = ThreadPoolExecutor(max_workers=3)
        threads.submit(send_frames)
        threads.submit(send_updates)
        threads.submit(recv_all)
        threads.shutdown(wait=True)


def get_server_port(server: Server):
    return server.socket.getsockname()[1]
