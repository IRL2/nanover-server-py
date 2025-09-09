import random
from concurrent.futures import ThreadPoolExecutor
from typing import Callable

import msgpack
from websockets.sync.client import connect

from nanover.state.state_dictionary import StateDictionary
from nanover.utilities.change_buffers import DictionaryChange


class WebsocketClient:
    def __init__(self, uri: str):
        self.current_frame = {}
        self.state_dictionary = StateDictionary()
        self.pending_commands = {}

        self.websocket = connect(uri)

        def listen():
            for message in self.websocket:
                self.recv_message(msgpack.unpackb(message))

        self.threads = ThreadPoolExecutor(max_workers=1)
        self.threads.submit(listen)

    def close(self):
        self.websocket.close()
        self.threads.shutdown(False)

    def run_command(self, name: str, arguments: dict, callback: Callable[[dict], None]):
        id = random.randint(0, 2**15)
        request = {"name": name, "arguments": arguments, "id": id}
        message = {"command": [{"request": request}]}
        self.pending_commands[id] = callback
        self.send_message(message)

    def send_message(self, message: dict):
        self.websocket.send(msgpack.packb(message))

    def recv_message(self, message: dict):
        if "frame" in message:
            self.recv_frame(message["frame"])
        if "state" in message:
            self.recv_state(message["state"])
        if "command" in message:
            for command in message["command"]:
                self.recv_command(command)

    def recv_frame(self, message: dict):
        self.current_frame.update(message)

    def recv_state(self, message: dict):
        change = DictionaryChange(
            updates=message["updates"],
            removals=message["removals"],
        )
        self.state_dictionary.update_state(None, change)

    def recv_command(self, message: dict):
        id = message["request"]["id"]
        response = message["response"]
        callback = self.pending_commands.pop(id)
        callback(response)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
