# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing a server for running a :class: ImdService.
"""
from typing import Optional

from narupa.core import NarupaServer, get_requested_port_or_default, DEFAULT_SERVE_ADDRESS
from narupa.imd.imd_service import ImdService
from narupa.protocol.imd.imd_pb2_grpc import add_InteractiveMolecularDynamicsServicer_to_server

DEFAULT_PORT = 54322


class ImdServer(NarupaServer):
    """
    Class providing an IMD server for running a :class: ImdService.

    :param: address: URL or IP address at which to run the server.
    :param: port: Port at which to run the server.
    """
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
        :return: The IMD service.
        """
        return self._imd_service

    def setup_services(self):
        """
        Sets up a new IMD service and attaches it to the server.
        """
        super().setup_services()
        self._imd_service = ImdService()
        self._imd_service.add_to_server_method(self._imd_service, self.server)
