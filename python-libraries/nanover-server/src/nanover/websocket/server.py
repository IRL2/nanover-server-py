import concurrent
import errno
from concurrent.futures import ThreadPoolExecutor
from ssl import SSLContext
from typing import Any, Self

import msgpack
import numpy as np

from nanover.core import AppServer
from nanover.trajectory import FrameData
from nanover.utilities.change_buffers import DictionaryChange
from nanover.utilities.cli import CancellationToken
from websockets.sync.server import serve, ServerConnection, Server

from nanover.websocket.commands import CommandMessageHandler


class WebSocketServer:
    @classmethod
    def basic_server(
        cls,
        app_server: AppServer,
        *,
        port: int = 0,
        ssl: SSLContext | None = None,
        insecure=True,
    ) -> "WebSocketServer":
        """
        Create a server for the given AppServer, listening for connection over one or both of insecure (ws) and
        secure (wss) websockets.
        """
        server = cls(app_server)

        if ssl is not None:
            server.serve(port=port, service="wss", ssl=ssl)
            port = 0  # choose random port from now on
        if insecure:
            server.serve(port=port, service="ws")

        return server

    def __init__(self, app_server: AppServer):
        self.app_server = app_server
        self._cancellation = CancellationToken()
        self._threads = ThreadPoolExecutor(thread_name_prefix="WebSocketServer")
        self._servers: list[Server] = []

    def close(self) -> None:
        """
        Close all open connections, subscriptions, and stop listening for new connections.
        """
        self._cancellation.cancel()

        for server in self._servers:
            server.shutdown()

        self._threads.shutdown()

    def serve(
        self,
        *,
        host="0.0.0.0",
        port=0,
        ssl: SSLContext | None = None,
        service: str | None = None,
    ) -> Server:
        """
        Listen for websocket connections on the given host and port, using SSL if a context is provided, and registering
        under discovery if a service name is provided.
        """
        try:
            server = serve(self._handle_client, host, port, ssl=ssl)
        except IOError as e:
            if e.errno == errno.EADDRINUSE:
                raise IOError(f"Port {port} already in use.") from None
            raise

        self._servers.append(server)
        self._threads.submit(server.serve_forever)
        if service:
            self.app_server.add_service(service, _get_server_port(server))
        return server

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

        self.user_id: str | None = None

        self._command_handler = CommandMessageHandler(
            websocket, command_service=app_server
        )

        self._threads = ThreadPoolExecutor(thread_name_prefix="WebSocketClientHandler")

    def close(self):
        print("CLOSING")
        self.cancellation.cancel()
        self._threads.shutdown()

    @property
    def frame_publisher(self):
        return self.app_server.frame_publisher

    @property
    def state_dictionary(self):
        return self.app_server.state_dictionary

    def send_frame(self, frame: FrameData):
        self.send_message({"frame": frame.pack_to_dict()})

    def send_state_update(self, change: DictionaryChange):
        self.send_message({"state": change.to_dict()})

    def send_message(self, message):
        self.websocket.send(msgpack.packb(message, default=_fallback_encoder))

    def recv_message(self, message: dict):
        def handle_state_update(update):
            change = DictionaryChange.from_dict(update)
            if self.user_id is None:
                self.user_id = _find_user_id(change)
            self.state_dictionary.update_state(None, change)

        if "state" in message:
            handle_state_update(message["state"])
        if "command" in message:
            # this can block waiting for another client to satisfy a command request
            self._threads.submit(
                self._command_handler.handle_message, message["command"]
            )

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

        concurrent.futures.wait(
            [
                self._threads.submit(send_frames),
                self._threads.submit(send_updates),
                self._threads.submit(recv_all),
            ]
        )
        self._threads.shutdown(wait=True)

        # remove keys owned by this user id
        if self.user_id is not None:
            removals = _find_user_owned_keys(self.app_server, self.user_id)
            self.state_dictionary.update_state(
                None, DictionaryChange(removals=removals)
            )


def _find_user_id(change: DictionaryChange):
    avatars = {key for key in change.updates if key.startswith("avatar.")}
    if avatars:
        return avatars.pop().replace("avatar.", "")
    return None


def _find_user_owned_keys(app_server: AppServer, user_id: str):
    """
    Return a set of state keys that are annotated as owned by the given user id.
    """
    keys = set()

    # user id as suffix
    with app_server.state_dictionary.lock_content() as state:
        keys |= {key for key in state if f".{user_id}" in key}

    # user id in interaction owner
    keys |= {
        key
        for key, interaction in app_server.imd.active_interactions.items()
        if interaction.properties.get("owner.id", None) == user_id
    }

    return keys


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
