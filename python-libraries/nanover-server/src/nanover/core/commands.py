from concurrent.futures import Future

from typing import Any, Protocol
from .app_server import CommandService as CommandServiceProtocol
from .types import CommandHandler, CommandRegistration


class SendMessage(Protocol):
    def __call__(self, message: Any) -> None: ...


class CommandMessageHandler:
    """
    Handle the message protocol for remote commands by handling incoming messages and emitting outgoing messages.
    This is essentially (part of) an RPC implementation.
    """

    def __init__(
        self,
        send_message: SendMessage,
        *,
        command_service: CommandServiceProtocol | None = None,
    ):
        self._command_service = command_service or CommandService()
        self._send_message = send_message

        self._pending_requests: dict[int, Future] = {}
        self._next_command_id = 1

    def register_command(
        self,
        name: str,
        callback: CommandHandler,
        default_arguments: dict | None = None,
    ) -> None:
        """Register a local callback that can be invoked by a remote party."""
        self._command_service.register_command(name, callback, default_arguments)
        self._send_message(
            {"register": {"name": name, "arguments": default_arguments}},
        )

    def request_command(
        self,
        name: str,
        arguments: dict | None = None,
    ) -> Future:
        """Request the invocation of a remote callback."""
        future: Future = Future()

        id = self._next_command_id
        self._next_command_id += 1
        request = {"name": name, "arguments": arguments or {}, "id": id}
        self._pending_requests[id] = future
        self._send_message({"request": request})
        return future

    def handle_message(self, message):
        """
        Handle an incoming command message, which may be a request to register a new command, a request to invoke a
        command, etc.

        Note: for now this will block when fulfilling a command registered remotely.
        """

        def handle_message(message: dict):
            if "register" in message:
                self._handle_command_registration(message["register"])
            elif "exception" in message:
                self._get_pending_future(message["request"]).set_exception(
                    RuntimeError(message["exception"])
                )
            elif "response" in message:
                self._get_pending_future(message["request"]).set_result(
                    message["response"]
                )
            else:
                self._handle_request(message["request"])

        # old format
        if isinstance(message, list):
            for submessage in message:
                handle_message(submessage)
        else:
            handle_message(message)

    def _get_pending_future(self, request):
        return self._pending_requests.pop(request["id"])

    def _handle_request(self, request):
        name, arguments = request.get("name"), request.get("arguments", {})

        def handle_future(future: Future):
            try:
                self._send_message({"request": request, "response": future.result()})
            except Exception as e:
                self._send_message({"request": request, "exception": str(e)})

        future = self._command_service.run_command(name, arguments)
        future.add_done_callback(handle_future)

    def _handle_command_registration(self, register):
        name, default_arguments = register.get("name"), register.get("arguments", {})

        def handle_call(**arguments):
            return self.request_command(name, arguments)

        self._command_service.register_command(name, handle_call, default_arguments)


class CommandService:
    """
    Implementation of the Command service, enabling services to register arbitrary commands
    which are run as callbacks.
    """

    def __init__(self, add_list_command=True):
        super().__init__()
        self.name: str = "command"
        self._commands = {}
        self._id = "service"

        def list_commands():
            return {
                "list": {
                    name: registration.arguments
                    for name, registration in self.commands.items()
                }
            }

        if add_list_command:
            self.register_command("commands/list", list_commands)

    @property
    def commands(self) -> dict[str, CommandRegistration]:
        """
        Gets a copy of the commands that have been registered, as :class:`CommandRegistration`,
        including their names, default arguments and registered callback.

        :return: A copy of the dictionary of commands that have been registered.
        """
        return self._commands.get_all()

    def register_command(
        self,
        name: str,
        callback: CommandHandler,
        default_arguments: dict | None = None,
    ):
        """
        Registers a command with this service

        :param name: Name of the command to register
        :param callback: Method to be called whenever the given command name is run by a client.
        :param default_arguments: A dictionary of the arguments of the callback and their default values.

        :raises ValueError: Raised when a command with the same name already exists.
        """
        if default_arguments is None:
            default_arguments = {}
        self._commands[name] = CommandRegistration(
            name=name,
            arguments=default_arguments,
            handler=callback,
        )

    def unregister_command(self, name):
        """
        Deletes a command from this service.

        :param name: Name of the command to delete
        """
        try:
            del self._commands[name]
        except KeyError:
            raise KeyError(f"Command {name} does not exist")

    def run_command(self, name: str, arguments: dict) -> Future:
        future: Future = Future()

        try:
            command = self._commands[name]

            if command is None:
                raise KeyError(f"Unknown command: {name}")

            result = command.run(arguments)
            if isinstance(result, Future):
                return result
        except Exception as e:
            future.set_exception(e)
        else:
            future.set_result(result)

        return future
