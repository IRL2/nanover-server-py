# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Reference multiplayer client implementation.

"""
import uuid
from typing import Optional
from narupa.core import NarupaClient


class MultiplayerClient(NarupaClient):
    """
    Represents a client to the multiplayer server.

    """
    _player_id: Optional[object] = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._player_id = str(uuid.uuid4())

    @property
    def player_id(self):
        """
        The player ID assigned to this client after joining multiplayer.

        :return: The player ID.
        """
        return self._player_id
