# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing an implementation of a multiplayer service,.
"""
import logging
from typing import Callable

from narupa.protocol.multiplayer.multiplayer_pb2 import (
    CreatePlayerRequest,
    CreatePlayerResponse,
)
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import (
    MultiplayerServicer,
    add_MultiplayerServicer_to_server,
)

MULTIPLAYER_SERVICE_NAME = "multiplayer"


class MultiplayerService(MultiplayerServicer):
    """
    Implementation of the Multiplayer service.
    """

    def __init__(self):
        super().__init__()
        self.name: str = MULTIPLAYER_SERVICE_NAME
        self.add_to_server_method: Callable = add_MultiplayerServicer_to_server

        self.players = {}
        self.logger = logging.getLogger(__name__)

    def CreatePlayer(self,
                     request: CreatePlayerRequest,
                     context) -> CreatePlayerResponse:
        """
        Create a new unique player and return their id.
        """
        player_id = self.generate_player_id()
        self.players[player_id] = request
        return CreatePlayerResponse(player_id=player_id)

    def generate_player_id(self):
        """
        Generates a new player ID.

        :return: A unique player ID.
        """
        return str(len(self.players) + 1)

    def close(self):
        pass
