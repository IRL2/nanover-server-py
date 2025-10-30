import time
from typing import Callable, Any, Protocol

from nanover.core.app_server import CommandService
from nanover.core.commands import CommandService as CommandServiceConcrete
from nanover.core.commands import CommandHandler


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
        command_service: CommandService | None = None
    ):
        self._command_service = command_service or CommandServiceConcrete()
        self._send_message = send_message

        self._pending_requests: dict[int, Callable[..., Any]] = {}
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
        callback: Callable[[dict], None] | None = None,
    ) -> None:
        """Request the invocation of a remote callback."""
        id = self._next_command_id
        self._next_command_id += 1
        request = {"name": name, "arguments": arguments or {}, "id": id}
        self._pending_requests[id] = callback or (lambda _: ...)
        self._send_message({"request": request})

    def handle_message(self, message):
        def handle_message(message):
            if "register" in message:
                self._handle_command_registration(message["register"])
            elif "exception" in message:
                self._resolve_pending_request(
                    message["request"], RuntimeError(message["exception"])
                )
            elif "response" in message:
                self._resolve_pending_request(message["request"], message["response"])
            else:
                self._handle_request(message["request"])

        # old format
        if isinstance(message, list):
            for submessage in message:
                handle_message(submessage)
        else:
            handle_message(message)

    def _resolve_pending_request(self, request, result):
        id = request["id"]
        callback = self._pending_requests.pop(id)
        callback(result)

    def _handle_request(self, request):
        name, arguments = request.get("name"), request.get("arguments", {})
        try:
            result = self._command_service.run_command(name, arguments)
            response = {"request": request, "response": result}
        except Exception as e:
            response = {"request": request, "exception": str(e)}
        self._send_message(response)

    def _handle_command_registration(self, register):
        name, arguments = register.get("name"), register.get("arguments", {})

        _UNRECEIVED = object()

        def handle_call_blocking(**arguments):
            returns = _UNRECEIVED

            def receive(results):
                nonlocal returns
                returns = results

            self.request_command(name, arguments, receive)

            while returns is _UNRECEIVED:
                time.sleep(0.1)

            if isinstance(returns, Exception):
                raise returns

            return returns

        self._command_service.register_command(name, handle_call_blocking, arguments)
