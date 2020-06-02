# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing an implementation of a multiplayer service,.
"""
import logging
from typing import Callable

MULTIPLAYER_SERVICE_NAME = "multiplayer"


class MultiplayerService:
    """
    Implementation of the Multiplayer service.
    """

    def __init__(self):
        super().__init__()
        self.name: str = MULTIPLAYER_SERVICE_NAME
        self.add_to_server_method: Callable = null

        self.players = {}
        self.logger = logging.getLogger(__name__)

    def generate_player_id(self):
        """
        Generates a new player ID.

        :return: A unique player ID.
        """
        return str(len(self.players) + 1)

    def close(self):
        pass


def null(*args, **kwargs):
    pass
