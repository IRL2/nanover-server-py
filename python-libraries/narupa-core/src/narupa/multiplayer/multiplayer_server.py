# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing implemenation of a muliplayer server.
"""
from typing import Optional

import grpc
from narupa.core import GrpcServer, get_requested_port_or_default, DEFAULT_SERVE_ADDRESS
from narupa.multiplayer.multiplayer_service import MultiplayerService
from narupa.protocol.multiplayer import multiplayer_pb2_grpc as multiplayer_proto_grpc

DEFAULT_PORT = 54323


class MultiplayerServer(GrpcServer):
    """
    Server providing multiplayer synchronisation.

    :param address: The IP or web address to run the server on.
    :param port: The port to run the server on.
    :param send_self: Whether to subscribe to own publications (e.g. receive own avatar).

    """

    multiplayer_service: MultiplayerService

    def __init__(self, *, address: Optional[str] = None,
                 port: Optional[int] = None):
        if address is None:
            address = DEFAULT_SERVE_ADDRESS
        port = get_requested_port_or_default(port, DEFAULT_PORT)
        self.address = address
        super().__init__(address=address, port=port)

    def setup_services(self):
        super().setup_services()
        self.multiplayer_service = MultiplayerService()
        multiplayer_proto_grpc.add_MultiplayerServicer_to_server(
            self.multiplayer_service, self.server)
