# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing an implementation of the :class:`CommandServicer`.

"""
from collections import namedtuple
from typing import Dict, Callable, Optional

import grpc
from google.protobuf.struct_pb2 import Struct

from narupa.core.command_info import CommandInfo
from narupa.core.protobuf_utilities import dict_to_struct
from narupa.core.key_lockable_map import KeyLockableMap
from narupa.protocol.command import CommandServicer, CommandMessage, CommandReply, GetCommandsReply

CommandRegistration = namedtuple('Command', ['info', 'callback'])


class CommandService(CommandServicer):
    """
    Implementation of the Command service, enabling services to register arbitrary commands
    which are run as callbacks.
    """

    def __init__(self):
        super().__init__()
        self._commands = KeyLockableMap()
        self._id = "service"

    @property
    def commands(self) -> Dict[str, CommandRegistration]:
        """
        Gets a copy of the commands that have been registered, as :class:`CommandRegistration`,
        including their names, default arguments and registered callback.

        :return: A copy of the dictionary of commands that have been registered.
        """
        return self._commands.get_all()

    def register_command(self, name: str, callback: Callable[[Dict[str, object]], Optional[Dict[str, object]]],
                         default_arguments: Optional[Dict[str, object]] = None):
        """
        Registers a command with this service

        :param name: Name of the command to register
        :param callback: Method to be called whenever the given command name is run by a client.
        :param default_arguments: A dictionary of the arguments of the callback and their default values.

        :raises ValueError: Raised when a command with the same name already exists.
        """
        if default_arguments is None:
            default_arguments = {}
        try:
            self._commands.set_no_replace(name, CommandRegistration(CommandInfo(name, **default_arguments), callback))
        except KeyError:
            raise ValueError(f"Command with name {name} has already been registered.")

    def unregister_command(self, name):
        """
        Deletes a command from this service.

        :param name: Name of the command to delete
        """
        if self._commands.get(name) is None:
            raise KeyError(f"Command {name} does not exist")
        self._commands.delete(self._id, name)

    def GetCommands(self, request, context):
        commands_copy = self.commands
        commands = [command.info.raw for command in
                    commands_copy.values()]
        return GetCommandsReply(commands=commands)

    def RunCommand(self, request, context):
        name = request.name
        command = self._commands.get(name)
        if command is None:
            context.set_code(grpc.StatusCode.INVALID_ARGUMENT)
            message = f'Unknown command: {command}'
            context.set_details(message)
            return
        args = command.info.arguments
        args.update(request.arguments)
        results = command.callback(**args)
        if results is not None:
            try:
                result_struct = dict_to_struct(results)
            except ValueError:
                context.set_code(grpc.StatusCode.INVALID_ARGUMENT)
                message = f'Command ({name}) generated results that cannot be serialised: {results}'
                context.set_details(message)
                return
            return CommandReply(result=result_struct)
        else:
            return CommandReply()
