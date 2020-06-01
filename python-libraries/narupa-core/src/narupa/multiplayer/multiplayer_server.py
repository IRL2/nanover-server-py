# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing an implementation of a multiplayer server.
"""
from typing import Optional

from narupa.core import DEFAULT_SERVE_ADDRESS
from narupa.core.grpc_server import (
    DEFAULT_MAX_WORKERS,
    get_requested_port_or_default,
)
from narupa.core import NarupaServer

DEFAULT_PORT = 54323


CREATE_ID_KEY = 'multiplayer/create_id'


class MultiplayerServer(NarupaServer):
    """
    Server providing multiplayer synchronisation.

    :param address: The IP or web address to run the server on.
    :param port: The port to run the server on.
    """

    def __init__(
            self,
            *,
            address: Optional[str] = None,
            port: Optional[int] = None,
            max_workers: int = DEFAULT_MAX_WORKERS,
    ):
        address = address or DEFAULT_SERVE_ADDRESS
        port = get_requested_port_or_default(port, DEFAULT_PORT)
        super().__init__(address=address, port=port, max_workers=max_workers)

    def setup_services(self):
        super().setup_services()
        next_id = 0

        def create_player_id():
            nonlocal next_id
            id = next_id
            next_id += 1
            return {'id': id}

        self.register_command(CREATE_ID_KEY, create_player_id)

