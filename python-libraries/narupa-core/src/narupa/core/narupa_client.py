# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Optional, Collection, Dict, List

from google.protobuf.struct_pb2 import Struct

from narupa.core import GrpcClient
from narupa.core.command_info import CommandInfo, struct_to_dict
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

    def run_command(self, name: str, **arguments) -> Dict[str, object]:
        """
        Runs a command on the command server.

        :param name: Name of command to run.
        :param arguments: Arguments to provide to command.

        :returns Dictionary of results, which may be empty.
        """
        arguments_struct = Struct()
        try:
            arguments_struct.update(arguments)
        except ValueError:
            raise ValueError("Unable to construct serialise arguments into a protobuf struct. "
                             "Only value types such as numbers, strings, booleans, and collections of those types"
                             "can be serialised.")

        message = CommandMessage(name=name, arguments=arguments_struct)
        result_message = self._command_stub.RunCommand(message)
        return struct_to_dict(result_message.result)

    def get_commands(self) -> List[CommandMessage]:
        """
        Get a list of all the commands on the command server. Blocks until the list of commands of
        available commands is received.

        :return: A list of all the commands on the command server.
        """
        return [CommandInfo.from_proto(raw) for raw in self._command_stub.GetCommands(GetCommandsRequest()).commands]