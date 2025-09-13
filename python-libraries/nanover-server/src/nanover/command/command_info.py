from typing import Dict

import msgpack

from nanover.utilities.protobuf_utilities import Serializable

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
        self._name = name
        self._arguments = arguments

    @property
    def name(self):
        return self._name

    @property
    def arguments(self):
        """
        Gets a copy of the default arguments this command accepts, as
        a dictionary.

        :return: Dictionary of default arguments.
        """
        return msgpack.unpackb(msgpack.packb(self._arguments))

    def __str__(self):
        args = self.arguments
        args_str = ", ".join(
            ["=".join([str(name), str(value)]) for name, value in args.items()]
        )
        return f"'{self.name}':({args_str})"
