# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Dict

from narupa.core.protobuf_utilities import dict_to_struct, struct_to_dict
from narupa.protocol.command import CommandMessage


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
    def arguments(self) -> Dict[str, object]:
        """
        Gets a copy of the default arguments this command accepts, as
        a dictionary.

        :return: Dictionary of default arguments.
        """
        return struct_to_dict(self.raw.arguments)
