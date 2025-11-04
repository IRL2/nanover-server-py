from concurrent.futures import Future

from nanover.core.app_server import CommandService
from nanover.core.types import CommandHandler, CommandRegistration
from nanover.websocket.client.base_client import WebsocketClient


class CommandClient(WebsocketClient, CommandService):
    """
    Mixin of methods for implementing CommandService
    """
    @property
    def commands(self) -> dict[str, CommandRegistration]:
        return self.run_command_blocking("commands/list")["list"]

    def run_command_blocking(self, name: str, **arguments):
        return self.run_command(name, arguments).result()

    def run_command(
        self,
        name: str,
        arguments: dict | None = None,
    ) -> Future:
        return self._command_handler.request_command(name, arguments)

    def register_command(
        self,
        name: str,
        callback: CommandHandler,
        default_arguments: dict | None = None,
    ) -> None:
        self._command_handler.register_command(name, callback, default_arguments)

    def unregister_command(self, name: str) -> None:
        raise NotImplementedError()
