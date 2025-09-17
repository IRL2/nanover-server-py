from typing import Any, Protocol

from nanover.app.commands import CommandHandler, CommandRegistration
from nanover.essd import ServiceHub
from nanover.imd import ImdStateWrapper
from nanover.state.state_dictionary import StateDictionary
from nanover.trajectory import FramePublisher
from nanover.utilities.change_buffers import DictionaryChange


class Closeable(Protocol):
    def close(self) -> None: ...


class StateService(Closeable, Protocol):
    def lock_state(self): ...

    def copy_state(self) -> dict[str, Any]: ...

    def update_state(self, access_token: Any, change: DictionaryChange): ...

    def clear_locks(self): ...

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
    ): ...

    def unregister_command(self, name: str): ...


class ImdService(Closeable, Protocol):
    @property
    def frame_publisher(self) -> FramePublisher: ...

    @property
    def imd(self) -> ImdStateWrapper: ...


class DiscoveryService(Closeable, Protocol):
    def add_service(self, name: str, port: int): ...

    @property
    def service_hub(self) -> ServiceHub: ...


class AppServer(CommandService, StateService, ImdService, DiscoveryService, Protocol):
    @property
    def name(self) -> str: ...
