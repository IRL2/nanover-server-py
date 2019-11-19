# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Optional, Dict

from google.protobuf.json_format import MessageToDict
from google.protobuf.struct_pb2 import Struct

from narupa.protocol.command import CommandMessage


def dict_to_struct(dictionary: Dict[str, object]) -> Struct:
    struct = Struct()
    try:
        struct.update(dictionary)
    except (ValueError, TypeError):
        raise ValueError("Unable to construct serialise dictionary into a protobuf struct. "
                         "Only value types such as numbers, strings, booleans, and collections of those types"
                         "can be serialised.")
    return struct


def struct_to_dict(struct: Struct) -> Dict[str, object]:
    return MessageToDict(struct)


class CommandInfo:
    """
    A wrapper around an underlying protobuf :class:`CommandMessage`,
    providing information about a given command.

    :param raw: The raw :class:`CommandMessage` this command is based on.
    """

    def __init__(self, name, **arguments):
        args_struct = dict_to_struct(arguments)
        self.raw = CommandMessage(name=name, arguments=args_struct)

    @classmethod
    def from_proto(cls, raw):
        args_raw = raw.arguments
        arguments = struct_to_dict(args_raw)
        return cls(raw.name, **arguments)

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
