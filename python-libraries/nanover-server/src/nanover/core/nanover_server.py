from typing import Callable, Optional, Dict, ContextManager, Set, Union

from nanover.command import CommandService
from nanover.command.command_service import CommandRegistration, CommandHandler
from nanover.core import GrpcServer
from nanover.state.state_service import StateService
from nanover.utilities.change_buffers import (
    DictionaryChange,
)

CommandCallable = Union[
    Callable[[Dict], Optional[Dict]], Callable[[], None], Callable[[Dict], None]
]


class NanoverServer(GrpcServer):
    """
    A base for NanoVer gRPC servers. Sets up a gRPC server, and automatically
    attaches a :class:`CommandService` and  :class:`StateService` enabling the running of arbitrary commands
    and synchronisation of state.
    """

    _command_service: CommandService
    _state_service: StateService

    def setup_services(self):
        """
        Sets up the services, including the :class:`CommandService`.
        """
        super().setup_services()
        self._setup_command_service()
        self._setup_state_service()

    def close(self):
        self._state_service.close()
        super().close()

    @property
    def commands(self) -> Dict[str, CommandRegistration]:
        """
        Gets the commands available on this server.

        :return: The commands, consisting of their names, callback and default parameters.
        """
        return self._command_service.commands

    def register_command(
        self,
        name: str,
        callback: CommandHandler,
        default_arguments: Optional[Dict] = None,
    ):
        """
        Registers a command with the :class:`CommandService` running on this server.

        :param name: Name of the command to register
        :param callback: Method to be called whenever the given command name is run by a client.
        :param default_arguments: A description of the arguments of the callback and their default values.

        :raises ValueError: Raised when a command with the same name already exists.
        """
        self._command_service.register_command(name, callback, default_arguments)

    def unregister_command(self, name):
        """
        Deletes a command from this service.

        :param name: Name of the command to delete
        """
        self._command_service.unregister_command(name)

    def lock_state(self) -> ContextManager[Dict[str, object]]:
        """
        Context manager for reading the current state while delaying any changes
        to it.
        """
        return self._state_service.lock_state()

    def copy_state(self) -> Dict[str, object]:
        """
        Return a deep copy of the current state.
        """
        return self._state_service.copy_state()

    def update_state(self, access_token: object, change: DictionaryChange):
        """
        Attempts an atomic update of the shared key/value store. If any key
        cannot be updates, no change will be made.
        """
        self._state_service.update_state(access_token, change)

    def update_locks(
        self,
        access_token: object = None,
        acquire: Optional[Dict[str, float]] = None,
        release: Optional[Set[str]] = None,
    ):
        """
        Attempts to acquire and release locks on keys in the shared key/value
        store. If any of the locks cannot be acquired, none of them will be.
        """
        if acquire is None:
            acquire = {}
        if release is None:
            release = set()
        self._state_service.update_locks(access_token, acquire, release)

    def clear_locks(self):
        """
        Forces the release all locks on all keys in the shared key/value store.
        """
        self._state_service.clear_locks()

    def _setup_command_service(self):
        self._command_service = CommandService()
        self.add_service(self._command_service)

    def _setup_state_service(self):
        self._state_service = StateService()
        self.add_service(self._state_service)
