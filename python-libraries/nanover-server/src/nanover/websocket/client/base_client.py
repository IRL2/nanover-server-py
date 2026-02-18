from concurrent.futures import ThreadPoolExecutor

import msgpack
from websockets.sync.client import connect, ClientConnection

from nanover.core.commands import CommandMessageHandler
from nanover.utilities.state_dictionary import StateDictionary
from nanover.utilities.change_buffers import DictionaryChange
from nanover.trajectory import FrameData
from nanover.websocket.server import MAX_MESSAGE_SIZE


class WebsocketClient:
    @classmethod
    def from_url(cls, url: str):
        connection = connect(url, max_size=MAX_MESSAGE_SIZE)
        client = cls(connection)
        return client

    def __init__(self, connection: ClientConnection):
        self._connection = connection

        def send_command(message):
            connection.send(msgpack.packb({"command": message}))

        self._state_dictionary = StateDictionary()
        self._command_handler = CommandMessageHandler(send_command)
        self._current_frame = FrameData()

        def listen():
            for bytes in self._connection:
                message = msgpack.unpackb(bytes)
                try:
                    self.recv_message(message)
                except Exception as e:
                    print(f"RECV FAILED ({set(message.keys())})", e)

        self.threads = ThreadPoolExecutor(thread_name_prefix="WebSocketClient")
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

    def send_message(self, message: dict):
        self._connection.send(msgpack.packb(message))

    def recv_message(self, message: dict):
        if "frame" in message:
            self.recv_frame(message["frame"])
        if "state" in message:
            self.recv_state(message["state"])
        if "command" in message:
            self._command_handler.handle_message(message["command"])

    def recv_frame(self, message: dict):
        self._current_frame.update(FrameData.unpack_from_dict(message))

    def recv_state(self, message: dict):
        change = DictionaryChange(
            updates=message["updates"],
            removals=message["removals"],
        )
        self._state_dictionary.update_state(None, change)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


_UNRECEIVED = object()
