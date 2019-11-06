# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Optional, Collection

from google.protobuf.struct_pb2 import Struct

from narupa.core import GrpcClient
from narupa.protocol.command import CommandStub, CommandMessage, GetCommandsRequest


class NarupaClient(GrpcClient):
    """
    A base gRPC client for Narupa services. Automatically sets up a stub for the :class:`CommandServicer`,
    enabling the running of arbitrary commands.

    :param address: Address of server to connect to.
    :param port: Port of server to connect to.
    :param stub: gRPC stub to attach.

    """

    def __init__(self, *, address: str,
                 port: int, stub: Optional = None):
        super().__init__(address=address, port=port,
                         stub=stub)

        self._command_stub = CommandStub(self.channel)

    def run_command(self, name: str, arguments: Optional[Struct] = None):
        """
        Runs a command on the command server.

        :param name: Name of command to run.
        :param arguments: Arguments to provide to command.
        """
        message = CommandMessage(name=name, arguments=arguments)
        return self._command_stub.RunCommand(message)

    def get_commands(self) -> Collection[CommandMessage]:
        """
        Get a list of all the commands on the command server.

        :return: A list of all the commands on the command server.
        """
        for reply in self._command_stub.GetCommands(GetCommandsRequest()):
            yield reply

    def get_commands_blocking(self):
        """
        Get a list of all the commands on the command server. Blocks until all commands received.

        :return: A list of all the commands on the command server.
        """
        return [c for c in self.get_commands()]
