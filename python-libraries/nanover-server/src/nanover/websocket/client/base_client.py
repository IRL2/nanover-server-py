import time
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Any

import msgpack
from websockets.sync.client import connect, ClientConnection

from nanover.state.state_dictionary import StateDictionary
from nanover.trajectory.frame_data import FRAME_INDEX
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.convert import unpack_dict_frame


MAX_MESSAGE_SIZE = 128 * 1024 * 1024


class WebsocketClient:
    @classmethod
    def from_url(cls, url: str):
        connection = connect(url, max_size=MAX_MESSAGE_SIZE)
        client = cls(connection)
        return client

    def __init__(self, connection: ClientConnection):
        self._connection = connection

        self._state_dictionary = StateDictionary()
        self._pending_commands: dict[int, Callable[..., Any]] = {}
        self._current_frame: dict[str, Any] = {}

        self.next_command_id = 1

        def listen():
            for bytes in self._connection:
                message = msgpack.unpackb(bytes)
                try:
                    self.recv_message(message)
                except Exception as e:
                    print(f"RECV FAILED ({set(message.keys())})", e)

        self.threads = ThreadPoolExecutor(
            max_workers=1, thread_name_prefix="WebSocketClient"
        )
        self.threads.submit(listen)

    def close(self):
        self._connection.close()
        self.threads.shutdown(True)

    def update_state(self, change):
        self.send_message(
            {
                "state": {
                    "updates": change.updates,
                    "removals": list(change.removals),
                }
            }
        )

    def run_command(
        self,
        name: str,
        arguments: dict | None = None,
        callback: Callable[[dict], None] | None = None,
    ):
        id = self.next_command_id
        self.next_command_id += 1
        request = {"name": name, "arguments": arguments or {}, "id": id}
        message = {"command": [{"request": request}]}
        self._pending_commands[id] = callback or (lambda _: ...)
        self.send_message(message)

    def run_command_blocking(self, name: str, **arguments):
        returns = _UNRECEIVED

        def receive(results):
            nonlocal returns
            returns = results

        self.run_command(name, arguments, receive)

        while returns is _UNRECEIVED:
            time.sleep(0.1)

        if isinstance(returns, Exception):
            raise returns

        return returns

    def send_message(self, message: dict):
        self._connection.send(msgpack.packb(message))

    def recv_message(self, message: dict):
        if "frame" in message:
            self.recv_frame(message["frame"])
        if "state" in message:
            self.recv_state(message["state"])
        if "command" in message:
            for command in message["command"]:
                self.recv_command(command)

    def recv_frame(self, message: dict):
        frame = unpack_dict_frame(message)
        if frame.get(FRAME_INDEX, None) == 0:
            self._current_frame = {}
        self._current_frame.update(frame)

    def recv_state(self, message: dict):
        change = DictionaryChange(
            updates=message["updates"],
            removals=message["removals"],
        )
        self._state_dictionary.update_state(None, change)

    def recv_command(self, message: dict):
        id = message["request"]["id"]

        if "exception" in message:
            response = RuntimeError(message["exception"])
        else:
            response = message["response"]

        print("CMD", response, message)

        callback = self._pending_commands.pop(id)
        callback(response)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


_UNRECEIVED = object()
