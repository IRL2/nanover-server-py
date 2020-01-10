# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

from uuid import uuid4
from typing import Dict, Iterable, ContextManager
from narupa.core import GrpcClient
from narupa.core.change_buffers import DictionaryChange
from narupa.core.command_info import CommandInfo
from narupa.core.protobuf_utilities import dict_to_struct, struct_to_dict, \
    deep_copy_dict
from narupa.core.state_dictionary import StateDictionary
from narupa.core.state_service import state_update_to_dictionary_change, \
    dictionary_change_to_state_update, validate_dict_is_serializable

from narupa.protocol.command import (
    CommandStub, CommandMessage, GetCommandsRequest,
)
from narupa.protocol.state import (
    StateStub, SubscribeStateUpdatesRequest, StateUpdate, UpdateStateRequest,
    UpdateLocksRequest,
)


DEFAULT_STATE_UPDATE_INTERVAL = 1 / 30


class NarupaClient(GrpcClient):
    """
    A base gRPC client for Narupa services. Automatically sets up a stub for the :class:`CommandServicer`,
    enabling the running of arbitrary commands.

    :param address: Address of server to connect to.
    :param port: Port of server to connect to.

    """
    _command_stub: CommandStub
    _state_stub: StateStub

    def __init__(self, *, address: str, port: int):
        super().__init__(address=address, port=port)

        self._command_stub = CommandStub(self.channel)
        self._available_commands = {}
        self._state_stub = StateStub(self.channel)
        self._access_token = str(uuid4())
        self._state = StateDictionary()

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
        arguments_struct = dict_to_struct(arguments)

        message = CommandMessage(name=name, arguments=arguments_struct)
        result_message = self._command_stub.RunCommand(message)
        return struct_to_dict(result_message.result)

    def update_available_commands(self) -> Dict[str, CommandInfo]:
        """
        Gets all the commands on the command server, and updates this
        client's known commands.
        Blocks until the dictionary of available commands is received.

        :return: A dictionary of all the commands on the command server, keyed by name
        """
        command_responses = self._command_stub.GetCommands(GetCommandsRequest()).commands
        self._available_commands = {raw.name: CommandInfo.from_proto(raw) for raw in command_responses}
        return self._available_commands

    def lock_state(self) -> ContextManager[Dict[str, object]]:
        """
        Context manager for reading the current state while delaying any changes
        to it.
        """
        return self._state.lock_content()

    def copy_state(self) -> Dict[str, object]:
        """
        Return a deep copy of the current state.
        """
        with self.lock_state() as state:
            return deep_copy_dict(state)

    def subscribe_all_state_updates(self, interval=DEFAULT_STATE_UPDATE_INTERVAL):
        """
        Subscribe, in the background, updates to the shared state.

        :param interval: Minimum time (in seconds) between receiving new updates
            for any and all values.
        """
        def process_state_updates(update_stream: Iterable[StateUpdate]):
            for update in update_stream:
                change = state_update_to_dictionary_change(update)
                self._state.update_state(None, change)

        request = SubscribeStateUpdatesRequest(update_interval=interval)
        update_stream = self._state_stub.SubscribeStateUpdates(request)
        self.threads.submit(process_state_updates, update_stream)

    def attempt_update_state(self, change: DictionaryChange) -> bool:
        validate_dict_is_serializable(change.updates)
        request = UpdateStateRequest(
            access_token=self._access_token,
            update=dictionary_change_to_state_update(change),
        )
        response = self._state_stub.UpdateState(request)
        return response.success

    def attempt_update_locks(self, **lock_updates: Dict[str, float]) -> bool:
        request = UpdateLocksRequest(
            access_token=self._access_token,
            key_changes=dict_to_struct(lock_updates),
        )
        response = self._state_stub.UpdateLocks(request)
        return response.success


class NarupaStubClient(NarupaClient):
    """
    A base gRPC client for Narupa services. Automatically sets up a stub for the :class:`CommandServicer`,
    and attaches the provided stub to the underlying gRPC channel.

    :param address: Address of server to connect to.
    :param port: Port of server to connect to.
    :param stub: gRPC stub to attach.

    """

    def __init__(self, *, address: str,
                 port: int, stub):
        super().__init__(address=address, port=port)
        self.stub = stub(self.channel)
