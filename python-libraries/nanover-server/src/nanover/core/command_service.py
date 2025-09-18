"""
Module providing an implementation of the :class:`CommandServicer`.

"""

from typing import Callable

import grpc

from nanover.core.commands import CommandService as CommandServiceBase
from nanover.protocol.command import (
    CommandServicer,
    CommandReply,
    CommandMessage,
    GetCommandsReply,
    add_CommandServicer_to_server,
)
from nanover.utilities.protobuf_utilities import dict_to_struct, struct_to_dict


class CommandService(CommandServicer, CommandServiceBase):
    """
    Implementation of the Command service, enabling services to register arbitrary commands
    which are run as callbacks.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.name: str = "command"
        self.add_to_server_method: Callable = add_CommandServicer_to_server
        self._id = "service"

    def GetCommands(self, request, context) -> GetCommandsReply:
        """
        GRPC method to get all of the commands available on this service.

        :param request: :class:`GetCommandsRequest`
        :param context: GRPC context.
        :return: :class:`GetCommandsReply`, detailing all the available commands.
        """
        commands = [
            CommandMessage(
                name=command.name, arguments=dict_to_struct(command.arguments)
            )
            for command in self.commands.values()
        ]
        return GetCommandsReply(commands=commands)

    def RunCommand(self, request, context) -> CommandReply:
        """
        GRPC method to run a command.

        :param request: :class:`CommandMessage` detailing the command to run and any arguments.
        :param context: GRPC context.
        :return: :class:`CommandReply`, consisting of any results of the command.
        """
        try:
            results = self.run_command(request.name, struct_to_dict(request.arguments))
            if results is None:
                return CommandReply()
            try:
                return CommandReply(result=dict_to_struct(results))
            except ValueError:
                context.set_code(grpc.StatusCode.INVALID_ARGUMENT)
                message = f"Command ({request.name}) generated results that cannot be serialised: {results}"
                context.set_details(message)
        except KeyError as e:
            context.set_code(grpc.StatusCode.INVALID_ARGUMENT)
            context.set_details(str(e))
