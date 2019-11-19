# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Optional, Collection, Dict, List, Set

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
        self._available_commands = Dict[str, CommandInfo]

    @property
    def available_commands(self) -> Dict[str, CommandInfo]:
        """
        Returns a copy of the dictionary of commands available on the server,
        as determined by previously calling :fun:`update_available_commands`.

        :return: A dictionary of command information, keyed by the command names.
        """
        return dict(self._available_commands)

    def run_command(self, name: str, **arguments) -> Dict[str, object]:
        """
        Runs a command on the command server.

        :param name: Name of command to run.
        :param arguments: Arguments to provide to command.

        :return: Dictionary of results, which may be empty.
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

    def update_available_commands(self) -> Dict[str, CommandInfo]:
        """
        Get a set of all the commands on the command server, and updates this
        clients set of known commands.
        Blocks until the list of commands of available commands is received.

        :return: A set of all the commands on the command server.
        """

        command_responses = self._command_stub.GetCommands(GetCommandsRequest()).commands
        self._available_commands = {raw.name: CommandInfo.from_proto(raw) for raw in command_responses}
        return self._available_commands
