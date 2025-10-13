from concurrent.futures import ThreadPoolExecutor
from ssl import SSLContext

import msgpack
import numpy as np
import numpy.typing as npt

from nanover.app.types import AppServer
from nanover.trajectory import FrameData2
from nanover.utilities.change_buffers import DictionaryChange
from nanover.utilities.cli import CancellationToken
from websockets.sync.server import serve, ServerConnection, Server
from nanover.trajectory.convert import pack_dict_frame

from types import TracebackType


DEFAULT_NANOVER_PORT = 38801


def validate_port(port: int) -> bool:
    """
    Ensures that the given `port` is valid.
    To be valid, the port must not use a protected address and must be accessible.

    :param port: Port ID to validate.
    :returns: True if valid port id else False.
    """
    if not isinstance(port, int):
        return False
    return port == 0 or (1023 < port < 65536)


class WebSocketServer:
    @classmethod
    def basic_server(
        cls,
        app_server: AppServer,
        *,
        port: int = 0,
        ssl: SSLContext | None = None,
        insecure: bool = True,
    ):
        server = cls(app_server)

        # TODO handle requesting the same port + proper reporting to user of the given port
        # maybe use logging or at the least have some reporter property
        if insecure:
            server.create_ws_server(port=port)
        if ssl is not None:
            port = (
                port + 1 if port != 0 else 0
            )  # Attempt to use next consecutive port if available.
            server.create_ws_server(port=port, ssl=ssl)

        return server

    def __init__(self, app_server: AppServer):
        self.app_server = app_server
        self._cancellation = CancellationToken()
        self._threads = ThreadPoolExecutor(
            max_workers=2, thread_name_prefix="WebSocketServer"
        )
        self._ws_server: Server | None = None
        self._wss_server: Server | None = None

    def create_ws_server(self, *, port: int = 0, ssl: SSLContext | None = None) -> int:
        """
        Creates WebSocket and attaches to the current object.
        Will attempt to create a Websocket at the desired `port`, using a random number (1023 < x < 65536)
        if it is already in use. If an SSLContext `ssl` is provided, a SSL wrapped socket
        is created instead.

        :param port: Prefered port number to use when creating the new WebSocket. A value of 0 will use a random port instead.
        :param ssl: `SSLContext` to use when verifying connections to the new WebSocket.
        :return: Port number for the new WebSocket server.
        """
        if not validate_port(port):
            raise ValueError(
                f"Invalid {port=}, please specify a number with the range 1024 - 65535 or 0 for a randomly assigned port address."
            )
        if ssl is None:
            type_, target = "ws", "_ws_server"
        else:
            type_, target = "wss", "_wss_server"

        if self.__getattribute__(target) is None:
            try:
                self.__setattr__(
                    target, serve(self._handle_client, "0.0.0.0", port, ssl=ssl)
                )
            except OSError as e:
                # OSError 98 indicates the port is already in use, so instead get a random available one.
                if e.errno != 98:
                    raise
                self.__setattr__(
                    target, serve(self._handle_client, "0.0.0.0", 0, ssl=ssl)
                )
            self._threads.submit(self.__getattribute__(target).serve_forever)
            self.app_server.add_service(type_, self.__getattribute__(type_ + "_port"))

        return get_server_port(self.__getattribute__(target))

    @property
    def ws_port(self) -> int | None:
        return get_server_port(self._ws_server) if self._ws_server is not None else None

    @property
    def wss_port(self) -> int | None:
        return (
            get_server_port(self._wss_server) if self._wss_server is not None else None
        )

    def close(self) -> None:
        self._cancellation.cancel()

        if self._wss_server is not None:
            self._wss_server.shutdown()
        if self._ws_server is not None:
            self._ws_server.shutdown()

        self._threads.shutdown()

    def __enter__(self) -> "WebSocketServer":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> bool:
        self.close()
        return False  # No attempt is made to handle exceptions so leave responsibility to caller.

    def _handle_client(self, websocket: ServerConnection) -> None:
        WebSocketClientHandler(self.app_server, websocket, self._cancellation).listen()


class WebSocketClientHandler:
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

    def send_frame(self, frame: FrameData2):
        self.send_message({"frame": pack_dict_frame(frame.frame_dict)})

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
        self.websocket.send(msgpack.packb(message, default=to_list))  # type: ignore

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


def get_server_port(server: Server) -> int:
    return server.socket.getsockname()[1]


def to_list(obj: npt.NDArray) -> list:
    """Converts `obj` to `list`, if of type `np.array` else raises `TypeError`."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Unknown type: {obj}")
