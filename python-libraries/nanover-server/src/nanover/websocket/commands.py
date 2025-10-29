from typing import Callable, Any

import msgpack
from websockets.sync.connection import Connection
from nanover.core.app_server import CommandService


class CommandMessageHandler:
    def __init__(self, command_service: CommandService, connection: Connection):
        self._command_service = command_service
        self._connection = connection

        self._pending_commands: dict[int, Callable[..., Any]] = {}
        self._next_command_id = 1

    def request_command(
        self,
        name: str,
        arguments: dict | None = None,
        callback: Callable[[dict], None] | None = None,
    ):
        id = self._next_command_id
        self._next_command_id += 1
        request = {"name": name, "arguments": arguments or {}, "id": id}
        self._pending_commands[id] = callback or (lambda _: ...)
        self.send_message({"request": request})

    def recv_message(self, message):
        def handle_message(message):
            if "exception" in message:
                self.handle_command_exception(message["request"], message["exception"])
            elif "response" in message:
                self.handle_command_response(message["request"], message["response"])
            else:
                self.handle_command_request(message["request"])

        # old format
        if isinstance(message, list):
            for submessage in message:
                handle_message(submessage)
        else:
            handle_message(message)

    def handle_command_exception(self, request, exception):
        id = request["id"]
        callback = self._pending_commands.pop(id)
        callback(RuntimeError(exception))

    def handle_command_response(self, request, response):
        id = request["id"]
        callback = self._pending_commands.pop(id)
        callback(response)

    def handle_command_request(self, request):
        name, arguments = request.get("name"), request.get("arguments", {})
        try:
            result = self._command_service.run_command(name, arguments)
            response = {"request": request, "response": result}
        except Exception as e:
            response = {"request": request, "exception": str(e)}
        self.send_message(response)

    def send_message(self, message):
        self._connection.send(msgpack.packb({"command": [message]}))
