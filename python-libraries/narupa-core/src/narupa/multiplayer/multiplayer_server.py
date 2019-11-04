# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing an implementation of a multiplayer server.
"""
from typing import Optional

from narupa.core import GrpcServer, DEFAULT_SERVE_ADDRESS
from narupa.core.grpc_server import DEFAULT_MAX_WORKERS, get_requested_port_or_default
from narupa.multiplayer.multiplayer_service import MultiplayerService
from narupa.protocol.multiplayer import multiplayer_pb2_grpc as multiplayer_proto_grpc

DEFAULT_PORT = 54323


class MultiplayerServer(GrpcServer):
    """
    Server providing multiplayer synchronisation.

    :param address: The IP or web address to run the server on.
    :param port: The port to run the server on.
    """
    _multiplayer_service: MultiplayerService

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
        self._multiplayer_service = MultiplayerService()
        multiplayer_proto_grpc.add_MultiplayerServicer_to_server(
            self._multiplayer_service, self.server)

    def close(self):
        super().close()
        self._multiplayer_service.close()
