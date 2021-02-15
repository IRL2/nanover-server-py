"""
A module for setting up typical Narupa clients, containing a client
that sets up a command service.
"""
# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

from typing import Dict, Iterable, ContextManager, Union, Any
from uuid import uuid4

import grpc
from narupa.command.command_info import CommandInfo
from narupa.core import GrpcClient
from narupa.protocol.command import (
    CommandStub, CommandMessage, GetCommandsRequest,
)
from narupa.protocol.state import (
    StateStub, SubscribeStateUpdatesRequest, StateUpdate, UpdateStateRequest,
    UpdateLocksRequest,
)
from narupa.state.state_dictionary import StateDictionary
from narupa.state.state_service import (
    state_update_to_dictionary_change, dictionary_change_to_state_update,
    validate_dict_is_serializable,
)
from narupa.utilities.change_buffers import DictionaryChange
from narupa.utilities.protobuf_utilities import (
    dict_to_struct, struct_to_dict, deep_copy_serializable_dict, Serializable,
)

DEFAULT_STATE_UPDATE_INTERVAL = 1 / 30


class NarupaClient(GrpcClient):
    """
    A base gRPC client for Narupa services. Automatically sets up a stub
    for the :class:`CommandServicer`, enabling the running of arbitrary commands.

    """
    _command_stub: CommandStub
    _available_commands: Dict[str, CommandInfo]
    _state_stub: StateStub
    _access_token: str
    _state: StateDictionary

    def __init__(self, *, channel: grpc.Channel, make_channel_owner: bool = False):
        super().__init__(channel=channel, make_channel_owner=make_channel_owner)
        self._setup_command_stub()
        self._setup_state_stub()

    def close(self):
        self._state.freeze()
        super().close()

    @property
    def available_commands(self) -> Dict[str, CommandInfo]:
        """
        Returns a copy of the dictionary of commands available on the server,
        as determined by previously calling :fun:`update_available_commands`.

        :return: A dictionary of command information, keyed by the command names.
        """
        return dict(self._available_commands)

    def run_command(self, name: str, **arguments) -> Dict[str, Serializable]:
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
        self._available_commands = {raw.name: CommandInfo.from_proto(raw)
                                    for raw in command_responses}
        return self._available_commands

    def lock_state(self) -> ContextManager[Dict[str, Serializable]]:
        """
        Context manager that locks and returns the state. Any attempted state
        updates are delayed until the context is exited.
        """
        return self._state.lock_content()

    def copy_state(self) -> Dict[str, Serializable]:
        """
        Return a deep copy of the current state.
        """
        with self.lock_state() as state:
            return deep_copy_serializable_dict(state)

    def subscribe_all_state_updates(self, interval=DEFAULT_STATE_UPDATE_INTERVAL):
        """
        Subscribe, in the background, to any updates made to the shared state.

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
        """
        Attempt to make a single atomic change to the shared state, blocking
        until a response is received.
        :param change: A single change to make to the shared state that will
            either be made in full, or ignored if some of the keys are locked
            by another user.
        :return: True if the server accepted our change, and False otherwise.
        """
        validate_dict_is_serializable(change.updates)
        request = UpdateStateRequest(
            access_token=self._access_token,
            update=dictionary_change_to_state_update(change),
        )
        response = self._state_stub.UpdateState(request)
        return response.success

    def attempt_update_locks(
            self,
            lock_updates: Dict[str, Union[float, None]]
    ) -> bool:
        """
        Attempt to acquire and/or free a number of locks on the shared state.
        :param lock_updates: A dictionary of keys to either a duration in
            seconds to attempt to acquire or renew a lock, or None to indicate
            the lock should be released if held.
        :return: True if the desired locks were acquired, and False otherwise.
        """
        request = UpdateLocksRequest(
            access_token=self._access_token,
            lock_keys=dict_to_struct(lock_updates),
        )
        response = self._state_stub.UpdateLocks(request)
        return response.success

    def _setup_command_stub(self):
        self._command_stub = CommandStub(self.channel)
        self._available_commands = {}

    def _setup_state_stub(self):
        self._state_stub = StateStub(self.channel)
        self._access_token = str(uuid4())
        self._state = StateDictionary()


class NarupaStubClient(NarupaClient):
    """
    A base gRPC client for Narupa services. Automatically sets up a stub
    for the :class:`CommandServicer`, and attaches the provided stub to
    the underlying gRPC channel.

    :param stub: gRPC stub to attach.

    """

    def __init__(self, *, channel: grpc.Channel, stub, make_channel_owner: bool = False):
        super().__init__(channel=channel, make_channel_owner=make_channel_owner)
        self.stub = stub(self.channel)
