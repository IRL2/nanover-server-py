# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing implemenation of a muliplayer server.
"""
from typing import Optional

import grpc
from narupa.core.grpc_server import  GrpcServer
from narupa.multiplayer.multiplayer_service import MultiplayerService
from narupa.protocol.multiplayer import multiplayer_pb2_grpc as multiplayer_proto_grpc

DEFAULT_ADDRESS = 'localhost'
DEFAULT_PORT = 7654


class MultiplayerServer(GrpcServer):
    """
    Server providing multiplayer synchronisation.

    :param address: The IP or web address to run the server on.
    :param port: The port to run the server on.
    :param send_self: Whether to subscribe to own publications (e.g. receive own avatar).

    """

    multiplayer_services: MultiplayerService

    def __init__(self, *, address: Optional[str] = None, port: Optional[int] =None, send_self=False):
        self.send_self = send_self
        if address is None:
            address = DEFAULT_ADDRESS
        if port is None:
            port = DEFAULT_PORT
        self.address = address
        self.port = port
        super().__init__(address=address, port=port)

    def setup_services(self):
        super().setup_services()
        self.multiplayer_services = MultiplayerService(send_self=self.send_self)
        multiplayer_proto_grpc.add_MultiplayerServicer_to_server(
            self.multiplayer_services, self.server)
