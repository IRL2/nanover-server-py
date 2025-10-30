from dataclasses import dataclass
from typing import Protocol
from nanover.utilities.key_lockable_map import KeyLockableMap


class CommandHandler(Protocol):
    def __call__(self, **kwargs) -> dict | None: ...


@dataclass(kw_only=True)
class CommandRegistration:
    name: str
    arguments: dict
    handler: CommandHandler

    def run(self, arguments: dict):
        args = {}
        args.update(self.arguments)
        args.update(arguments)
        return self.handler(**args)


class CommandService:
    """
    Implementation of the Command service, enabling services to register arbitrary commands
    which are run as callbacks.
    """

    def __init__(self, add_list_command=True):
        super().__init__()
        self.name: str = "command"
        self._commands = KeyLockableMap()
        self._id = "service"

        def list_commands():
            return {
                "list": {
                    name: registration.arguments
                    for name, registration in self.commands.items()
                }
            }

        if add_list_command:
            self.register_command("commands/list", list_commands)

    @property
    def commands(self) -> dict[str, CommandRegistration]:
        """
        Gets a copy of the commands that have been registered, as :class:`CommandRegistration`,
        including their names, default arguments and registered callback.

        :return: A copy of the dictionary of commands that have been registered.
        """
        return self._commands.get_all()

    def register_command(
        self,
        name: str,
        callback: CommandHandler,
        default_arguments: dict | None = None,
    ):
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
            self._commands.set_no_replace(
                name,
                CommandRegistration(
                    name=name,
                    arguments=default_arguments,
                    handler=callback,
                ),
            )
        except KeyError:
            raise ValueError(f"Command with name {name} has already been registered.")

    def unregister_command(self, name):
        """
        Deletes a command from this service.

        :param name: Name of the command to delete
        """
        try:
            self._commands.delete(self._id, name)
        except KeyError:
            raise KeyError(f"Command {name} does not exist")

    def run_command(self, name: str, arguments: dict):
        command = self._commands.get(name)
        if command is None:
            raise KeyError(f"Unknown command: {name}")
        return command.run(arguments)
