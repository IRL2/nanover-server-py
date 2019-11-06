# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Callable, Optional, Dict

from google.protobuf.struct_pb2 import Struct

from narupa.command import CommandService, Command
from narupa.core import GrpcServer
from narupa.core.grpc_server import DEFAULT_MAX_WORKERS
from narupa.protocol.command import add_CommandServicer_to_server


class NarupaServer(GrpcServer):
    """
    A base for Narupa gRPC servers. Sets up a gRPC server, and automatically
    attaches a :class:`CommandService`, enabling the running of arbitrary commands.
    """

    def __init__(self, *, address: str, port: int, max_workers=DEFAULT_MAX_WORKERS):
        super().__init__(address=address, port=port, max_workers=max_workers)

    def setup_services(self):
        """
        Sets up the :class:`CommandService`.
        """
        super().setup_services()
        self._command_service = CommandService()
        add_CommandServicer_to_server(self._command_service, self.server)

    @property
    def commands(self) -> Dict[str, Command]:
        """
        Gets the commands available on this server.

        :return: The commands, consisting of their names, callback and default parameters.
        """
        return self._command_service.commands

    def register_command(self, name: str, callback: Callable[[Struct], Optional[Struct]],
                         default_arguments: Optional[Struct] = None):
        """
        Registers a command with the :class:`CommandService` running on this server.

        :param name: Name of the command to register
        :param callback: Method to be called whenever the given command name is run by a client.
        :param default_arguments: A description of the arguments of the callback and their default values.

        :raises: ValueError: Raised when a command with the same name already exists.
        """
        self._command_service.register_command(name, callback, default_arguments)

    def delete_command(self, name):
        """
        Deletes a command from this service.

        :param name: Name of the command to delete
        """
        self._command_service.delete_command(name)
