from typing import Any, Protocol

from nanover.core.commands import CommandRegistration, CommandHandler
from nanover.essd import ServiceHub, DiscoveryServer
from nanover.imd import ImdStateWrapper
from nanover.utilities.state_dictionary import StateDictionary
from nanover.trajectory import FramePublisher
from nanover.utilities.change_buffers import DictionaryChange


class Closeable(Protocol):
    def close(self) -> None: ...


class StateService(Closeable, Protocol):
    def lock_state(self) -> None: ...

    def copy_state(self) -> dict[str, Any]: ...

    def update_state(self, access_token: Any, change: DictionaryChange) -> None: ...

    def clear_locks(self) -> None: ...

    @property
    def state_dictionary(self) -> StateDictionary: ...


class CommandService(Protocol):
    @property
    def commands(self) -> dict[str, CommandRegistration]: ...

    def run_command(self, name: str, arguments: dict[str, Any]) -> dict[str, Any]: ...

    def register_command(
        self,
        name: str,
        callback: CommandHandler,
        default_arguments: dict | None = None,
    ) -> None: ...

    def unregister_command(self, name: str) -> None: ...


class ImdService(Closeable, Protocol):
    @property
    def frame_publisher(self) -> FramePublisher: ...

    @property
    def imd(self) -> ImdStateWrapper: ...


class DiscoveryService(Closeable, Protocol):
    @property
    def discovery(self) -> DiscoveryServer | None: ...

    def add_service(self, name: str, port: int) -> None: ...

    @property
    def service_hub(self) -> ServiceHub: ...


class AppServerMinimal(CommandService, StateService, Protocol): ...


class AppServer(CommandService, StateService, ImdService, DiscoveryService, Protocol):
    @property
    def name(self) -> str: ...


def basic_info_string(app_server: AppServer):
    protocols = ", ".join(
        f"{protocol}://localhost:{port}"
        for protocol, port in app_server.service_hub.services.items()
    )

    discovery = app_server.discovery

    return f'Serving "{app_server.name}" ({protocols}), ' + (
        f"discoverable on all interfaces on port {discovery.port}"
        if discovery
        else "without discovery"
    )
