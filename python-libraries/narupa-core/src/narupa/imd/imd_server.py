# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

from typing import Optional

from narupa.core import GrpcServer, get_requested_port_or_default, DEFAULT_SERVE_ADDRESS
from narupa.imd.imd_service import ImdService
from narupa.protocol.imd.imd_pb2_grpc import add_InteractiveMolecularDynamicsServicer_to_server

DEFAULT_PORT = 54322


class ImdServer(GrpcServer):
    _imd_service: ImdService

    def __init__(self, *, address: Optional[str]=None, port: Optional[int]=None):
        if address is None:
            address = DEFAULT_SERVE_ADDRESS
        port = get_requested_port_or_default(port, DEFAULT_PORT)
        super().__init__(address=address, port=port)

    @property
    def service(self) -> ImdService:
        """
        Gets the IMD service implementation attached to this server.
        :return:
        """
        return self._imd_service

    def setup_services(self):
        super().setup_services()
        self._imd_service = ImdService()
        add_InteractiveMolecularDynamicsServicer_to_server(self._imd_service, self.server)
