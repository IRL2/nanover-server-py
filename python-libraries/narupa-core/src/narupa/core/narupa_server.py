# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Callable, Optional, Dict, ContextManager, Set

import grpc
from google.protobuf.struct_pb2 import Struct
from narupa.core.change_buffers import DictionaryChange

from narupa.core.command_service import CommandService, CommandRegistration
from narupa.core.state_service import StateService
from narupa.core import GrpcServer
from narupa.core.grpc_server import DEFAULT_MAX_WORKERS
from narupa.protocol.command import add_CommandServicer_to_server
from narupa.protocol.state import add_StateServicer_to_server


class NarupaServer(GrpcServer):
    """
    A base for Narupa gRPC servers. Sets up a gRPC server, and automatically
    attaches a :class:`CommandService`, enabling the running of arbitrary commands.
    """
    _command_service: CommandService
    _state_service: StateService

    def __init__(self, *, address: str, port: int, max_workers=DEFAULT_MAX_WORKERS):
        super().__init__(address=address, port=port, max_workers=max_workers)

    def setup_services(self):
        """
        Sets up the services, including the :class:`CommandService`.
        """
        super().setup_services()
        self._setup_command_service()
        self._setup_state_service()

    @property
    def commands(self) -> Dict[str, CommandRegistration]:
        """
        Gets the commands available on this server.

        :return: The commands, consisting of their names, callback and default parameters.
        """
        return self._command_service.commands

    def register_command(self, name: str, callback: Callable[[Dict], Optional[Dict]],
                         default_arguments: Optional[Dict] = None):
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
            acquire: Dict[str, float] = {},
            release: Set[str] = set(),
    ):
        """
        Attempts to acquire and release locks on keys in the shared key/value
        store. If any of the locks cannot be acquired, none of them will be.
        """
        self._state_service.update_locks(access_token, acquire, release)

    def _setup_command_service(self):
        self._command_service = CommandService()
        self.add_service(self._command_service, add_CommandServicer_to_server)


    def _setup_state_service(self):
        self._state_service = StateService()
        self.add_service(self._state_service, add_StateServicer_to_server)
