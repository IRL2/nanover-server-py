# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing an implementation of a client for :class:`CommandServicer`, primarily for testing.
"""
from typing import Optional, Collection

from google.protobuf.struct_pb2 import Struct

from narupa.command.command_server import DEFAULT_PORT
from narupa.core import GrpcClient, get_requested_port_or_default
from narupa.protocol.command import CommandStub, CommandMessage, GetCommandsRequest


class CommandClient(GrpcClient):
    """
    A gRPC client for the :class:`CommandServicer`.

    :param address: Address of server to connect to.
    :param port: Port of server to connect to.

    """

    def __init__(self, *, address: Optional[str] = None,
                 port: Optional[int] = None):
        port = get_requested_port_or_default(port, DEFAULT_PORT)
        super().__init__(address=address, port=port,
                         stub=CommandStub)

    def run_command(self, name: str, arguments: Optional[Struct] = None):
        """
        Runs a command on the command server.

        :param name: Name of command to run.
        :param arguments: Arguments to provide to command.
        """
        if arguments is None:
            arguments = Struct()
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
