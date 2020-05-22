# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Reference multiplayer client implementation.

"""

from typing import Callable, Sequence

import grpc
import narupa.protocol.multiplayer.multiplayer_pb2_grpc as mult_proto_grpc
from narupa.core import NarupaStubClient
from narupa.multiplayer.multiplayer_server import DEFAULT_PORT
from narupa.protocol.multiplayer.multiplayer_pb2 import CreatePlayerRequest
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import MultiplayerStub
from narupa.utilities.protobuf_utilities import struct_to_dict

UpdateCallback = Callable[[Sequence[str]], None]


def _end_upon_channel_close(function):
    """
    Wrapper function to streaming RPCs that gracefully handles any exceptions thrown when the channel has closed.
    :param function:
    :return:
    """

    def wrapped(self, *args, **kwargs):
        try:
            function(self, *args, **kwargs)
        except grpc.RpcError as e:
            if self.closed:
                return
            if e.args and 'Channel closed' in e.args[0]:
                return
            else:
                raise e

    return wrapped


class MultiplayerClient(NarupaStubClient):
    """
    Represents a client to the multiplayer server.

    """
    stub: mult_proto_grpc.MultiplayerStub
    DEFAULT_CONNECTION_PORT = DEFAULT_PORT

    def __init__(self, *,
                 channel: grpc.Channel,
                 make_channel_owner: bool = False):
        super().__init__(channel=channel,
                         stub=MultiplayerStub,
                         make_channel_owner=make_channel_owner)
        self._player_id = None

    def join_multiplayer(self, player_name):
        """
        Joins a multiplayer server

        :param player_name: The user name to use for multiplayer.
        :return: Player ID received from the multiplayer server.
        """
        if self.joined_multiplayer:
            return self._player_id
        response = self.stub.CreatePlayer(CreatePlayerRequest(player_name=player_name), timeout=5)
        self._player_id = response.player_id
        return self._player_id

    @property
    def player_id(self):
        """
        The player ID assigned to this client after joining multiplayer.

        :return: The player ID.
        """
        return self._player_id

    @property
    def joined_multiplayer(self) -> bool:
        """
        Indicates whether multiplayer joined, and the client has a valid player ID.

        :return: True if multiplayer has been joined, false otherwise.
        """
        return self.player_id is not None

    def _assert_has_player_id(self):
        if not self.joined_multiplayer:
            raise RuntimeError("Join multiplayer before attempting this operation.")
