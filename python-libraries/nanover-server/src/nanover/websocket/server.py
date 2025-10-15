from concurrent.futures import ThreadPoolExecutor
from ssl import SSLContext
from typing import Any, Self

import msgpack
import numpy as np

from nanover.app.types import AppServer
from nanover.trajectory import FrameData
from nanover.utilities.change_buffers import DictionaryChange
from nanover.utilities.cli import CancellationToken
from websockets.sync.server import serve, ServerConnection, Server


class WebSocketServer:
    @classmethod
    def basic_server(
        cls,
        app_server: AppServer,
        *,
        ssl: SSLContext | None = None,
        insecure=True,
    ) -> "WebSocketServer":
        """
        Create a server for the given AppServer, listening for connection over one or both of insecure (ws) and
        secure (wss) websockets.
        """
        server = cls(app_server)

        if insecure:
            server.serve_insecure()
        if ssl is not None:
            server.serve_secure(ssl=ssl)

        return server

    def __init__(self, app_server: AppServer):
        self.app_server = app_server
        self._cancellation = CancellationToken()
        self._threads = ThreadPoolExecutor(
            max_workers=2, thread_name_prefix="WebSocketServer"
        )
        self._ws_server: Server | None = None
        self._wss_server: Server | None = None

    def close(self) -> None:
        """
        Close all open connections, subscriptions, and stop listening for new connections.
        """
        self._cancellation.cancel()

        if self._wss_server is not None:
            self._wss_server.shutdown()
        if self._ws_server is not None:
            self._ws_server.shutdown()

        self._threads.shutdown()

    def serve_insecure(self, *, host="0.0.0.0", port=0) -> Server:
        """
        Listen for insecure websocket (ws) connections on the given host and port.
        """
        if self._ws_server is None:
            self._ws_server = serve(self._handle_client, host, port)
            self._threads.submit(self._ws_server.serve_forever)
            self.app_server.add_service("ws", _get_server_port(self._ws_server))
        return self._ws_server

    def serve_secure(self, *, host="0.0.0.0", port=0, ssl: SSLContext) -> Server:
        """
        Listen for secure websocket (wss) connections on the given host and port using the given SSL context.
        """
        if self._wss_server is None:
            self._wss_server = serve(self._handle_client, host, port, ssl=ssl)
            self._threads.submit(self._wss_server.serve_forever)
            self.app_server.add_service("wss", _get_server_port(self._wss_server))
        return self._wss_server

    def _handle_client(self, websocket: ServerConnection) -> None:
        _WebSocketClientHandler(self.app_server, websocket, self._cancellation).listen()

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class _WebSocketClientHandler:
    """
    Handles the connection of a single client to the server, subscribing to the server's app server and passing on
    state and frame updates, as well as handling incoming state updates and command requests.
    """

    def __init__(
        self,
        app_server: AppServer,
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
        self.send_message({"frame": frame.pack_to_dict()})

    def send_state_update(self, change: DictionaryChange):
        self.send_message({"state": change.to_dict()})

    def send_message(self, message):
        self.websocket.send(msgpack.packb(message, default=_fallback_encoder))

    def recv_message(self, message: dict):
        def handle_state_update(update):
            self.state_dictionary.update_state(None, DictionaryChange.from_dict(update))

        def handle_command_request(request):
            name, arguments = request.get("name"), request.get("arguments", {})
            try:
                result = self.run_command(name, arguments)
                response = {"request": request, "response": result}
            except Exception as e:
                response = {"request": request, "exception": str(e)}
            self.send_message({"command": [response]})

        if "state" in message:
            handle_state_update(message["state"])
        if "command" in message:
            # old format
            if isinstance(message["command"], list):
                for submessage in message["command"]:
                    handle_command_request(submessage["request"])
            else:
                handle_command_request(message["command"]["request"])

    def listen(self, frame_interval=1 / 30, state_interval=1 / 30):
        # TODO: error handling!!
        def send_frames():
            for frame in self.frame_publisher.subscribe_latest_frames(
                frame_interval=frame_interval,
                cancellation=self.cancellation,
            ):
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


def _get_server_port(server: Server) -> int:
    """
    Returns the concrete port used by a given websocket server.
    """
    return server.socket.getsockname()[1]


def _fallback_encoder(obj: Any) -> Any:
    """
    Converts, if possible, a type msgpack doesn't understand into a basic type it can encode.

    :param obj: object to be converted
    :return: simplified object
    """
    # encode numpy arrays as simple lists
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Unknown type: {obj}")
