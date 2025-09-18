from concurrent.futures import ThreadPoolExecutor
from ssl import SSLContext

import msgpack

from nanover.app import NanoverImdApplication
from nanover.trajectory.frame_data import FrameData, FRAME_INDEX
from nanover.utilities.change_buffers import DictionaryChange
from nanover.utilities.cli import CancellationToken
from websockets.sync.server import serve, ServerConnection, Server

from nanover.websocket.convert import pack_grpc_frame


class WebSocketServer:
    @classmethod
    def basic_server(
        cls,
        app_server: NanoverImdApplication,
        *,
        ssl: SSLContext | None = None,
        insecure=True,
    ):
        server = cls(app_server)

        if insecure:
            server.serve_insecure()
        if ssl is not None:
            server.serve_secure(ssl=ssl)

        return server

    def __init__(self, app_server: NanoverImdApplication):
        self.app_server = app_server
        self._cancellation = CancellationToken()
        self._threads = ThreadPoolExecutor(
            max_workers=2, thread_name_prefix="WebSocketServer"
        )
        self._ws_server: Server | None = None
        self._wss_server: Server | None = None

    def serve_insecure(self):
        if self._ws_server is None:
            self._ws_server = serve(self._handle_client, "0.0.0.0", 0)
            self._threads.submit(self._ws_server.serve_forever)
            self.app_server.add_service("ws", self.ws_port)

    def serve_secure(self, *, ssl: SSLContext):
        if self._wss_server is None:
            self._wss_server = serve(self._handle_client, "0.0.0.0", 0, ssl=ssl)
            self._threads.submit(self._wss_server.serve_forever)
            self.app_server.add_service("wss", self.wss_port)

    @property
    def ws_port(self):
        return get_server_port(self._ws_server) if self._ws_server is not None else None

    @property
    def wss_port(self):
        return (
            get_server_port(self._wss_server) if self._wss_server is not None else None
        )

    def close(self):
        self._cancellation.cancel()

        if self._wss_server is not None:
            self._wss_server.shutdown()
        if self._ws_server is not None:
            self._ws_server.shutdown()

        self._threads.shutdown()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _handle_client(self, websocket: ServerConnection):
        WebSocketClientHandler(self.app_server, websocket, self._cancellation).listen()


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
        cancellation.subscribe_cancellation(self.close)

    def close(self):
        self.cancellation.cancel()

    @property
    def frame_publisher(self):
        return self.app_server.frame_publisher

    @property
    def state_dictionary(self):
        return self.app_server.state_dictionary

    def run_command(self, name: str, arguments: dict | None = None):
        results = self.app_server.run_command(name, arguments or {})
        return results

    def send_frame(self, frame: FrameData):
        self.send_message({"frame": pack_grpc_frame(frame)})

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
            try:
                result = self.run_command(name, arguments)
                response = {"request": request, "response": result}
            except Exception as e:
                response = {"request": request, "exception": str(e)}
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
                frame = FrameData(response.frame)
                frame.values[FRAME_INDEX] = response.frame_index
                self.send_frame(frame)

        def send_updates():
            with self.state_dictionary.get_change_buffer() as change_buffer:
                self.cancellation.subscribe_cancellation(change_buffer.freeze)
                for change in change_buffer.subscribe_changes(interval=state_interval):
                    self.send_state_update(change)

        def recv_all():
            for data in self.websocket:
                self.recv_message(msgpack.unpackb(data))
            self.cancellation.cancel()

        threads = ThreadPoolExecutor(
            max_workers=3, thread_name_prefix="WebSocketClientHandler"
        )
        threads.submit(send_frames)
        threads.submit(send_updates)
        threads.submit(recv_all)
        threads.shutdown(wait=True)


def get_server_port(server: Server):
    return server.socket.getsockname()[1]
