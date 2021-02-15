# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Dict

from narupa.utilities.protobuf_utilities import (
    dict_to_struct, struct_to_dict, Serializable,
)
from narupa.protocol.command import CommandMessage

CommandArguments = Dict[str, Serializable]
CommandResult = Dict[str, Serializable]


class CommandInfo:
    """
    A wrapper around an underlying protobuf :class:`CommandMessage`,
    providing information about a given command.

    :param name: Name of the command.
    :param arguments: Dictionary of command arguments.
    """

    def __init__(self, name, **arguments):
        args_struct = dict_to_struct(arguments)
        self.raw = CommandMessage(name=name, arguments=args_struct)

    @classmethod
    def from_proto(cls, raw):
        instance = cls(raw.name)
        instance.raw.MergeFrom(raw)
        return instance

    @property
    def name(self) -> str:
        return self.raw.name

    @property
    def arguments(self) -> CommandArguments:
        """
        Gets a copy of the default arguments this command accepts, as
        a dictionary.

        :return: Dictionary of default arguments.
        """
        return struct_to_dict(self.raw.arguments)

    def __str__(self):
        args = self.arguments
        args_str = ', '.join(['='.join([str(name), str(value)]) for name, value in args.items()])
        return f'\'{self.name}\':({args_str})'
