import msgpack
from websockets.sync.connection import Connection

from nanover.core import AppServer


class CommandMessageHandler:
    def __init__(self, app_server: AppServer, connection: Connection):
        self._app_server = app_server
        self._connection = connection

    def recv_message(self, message):
        def handle_message(message):
            if "response" in message:
                self.handle_command_response(message["request"], message["response"])

            self.handle_command_request(message["request"])

        # old format
        if isinstance(message, list):
            for submessage in message:
                handle_message(submessage)
        else:
            handle_message(message)

    def handle_command_response(self, request, response):
        pass

    def handle_command_request(self, request):
        name, arguments = request.get("name"), request.get("arguments", {})
        try:
            result = self._app_server.run_command(name, arguments)
            response = {"request": request, "response": result}
        except Exception as e:
            response = {"request": request, "exception": str(e)}
        self.send_message(response)

    def send_message(self, message):
        self._connection.send(msgpack.packb({"command": [message]}))
