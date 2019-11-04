from collections import namedtuple
from typing import Dict, Callable, Optional

import grpc
from google.protobuf.struct_pb2 import Struct

from narupa.multiplayer.key_lockable_map import KeyLockableMap
from narupa.protocol.command import CommandServicer, CommandMessage, CommandReply

Command = namedtuple('Command', ['callback', 'default_args'])


class CommandService(CommandServicer):
    """
    Implementation of the Command service, enabling services to register arbitrary commands.
    """

    def __init__(self):
        super().__init__()
        self._commands = KeyLockableMap()
        self._id = "service"

    @property
    def commands(self) -> Dict[str, Command]:
        return self._commands.get_all()

    def register_command(self, name: str, callback: Callable[[Struct], Optional[Struct]], default_arguments: Struct):
        """
        Registers a command with this service

        :param name: Name of the command to register
        :param callback: Method to be called whenever this command is run by a client.

        :raises: ValueError: Raised when a command with the same name already exists.
        """
        if self._commands.get(name) is not None:
            raise ValueError(f"Command with name {name} has already been registered.")
        self._commands.set(self._id, name, Command(callback, default_arguments))

    def delete_command(self, name):
        """
        Deletes a command from this service.

        :param name: Name of the command to delete
        """
        self._commands.delete(self._id, name)

    def GetCommands(self, request, context):
        commands_copy = self.commands
        for name, command in commands_copy.items():
            command_message = CommandMessage(name=name, arguments=command.default_args)
            yield command_message

    def RunCommand(self, request, context):
        name = request.name
        command = self._commands.get(name)
        if command is None:
            context.set_code(grpc.StatusCode.INVALID_ARGUMENT)
            message = f'Unknown command: {command}'
            context.set_details(message)
            return
        results = command.callback(request.arguments)
        if type(results) is Struct:
            return CommandReply(result=results)
        else:
            return CommandReply()
