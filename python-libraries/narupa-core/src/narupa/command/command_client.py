"""
Implementation of a client for :class:`CommandServicer`, primarily for testing.
"""
from typing import Optional, Collection

from google.protobuf.internal.well_known_types import Struct

from narupa.command.command_server import DEFAULT_PORT
from narupa.core import GrpcClient, get_requested_port_or_default
from narupa.protocol.command import CommandStub, CommandMessage, GetCommandsRequest


class CommandClient(GrpcClient):

    def __init__(self, *, address: Optional[str] = None,
                 port: Optional[int] = None):
        port = get_requested_port_or_default(port, DEFAULT_PORT)
        super().__init__(address=address, port=port,
                         stub=CommandStub)

    def run_command(self, name: str, arguments: Struct):
        """
        Runs a command on the command server.

        :param name: Name of command to run.
        :param arguments: Arguments to provide to command.
        """
        message = CommandMessage(name=name, arguments=arguments)
        return self.stub.RunCommand(message)

    def get_commands(self) -> Collection[CommandMessage]:
        """
        Get a list of all the commands on the command server.

        :return: A list of all the commands on the command server.
        """
        for reply in self.stub.GetCommands(GetCommandsRequest()):
            yield reply

    def get_commands_blocking(self):
        """
        Get a list of all the commands on the command server. Blocks until all commands received.

        :return: A list of all the commands on the command server.
        """
        return [c for c in self.get_commands()]