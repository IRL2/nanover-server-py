# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing an implementation of a server for :class:`CommandServicer`, primarily for testing.
"""
from typing import Optional

from narupa.core import GrpcServer, DEFAULT_SERVE_ADDRESS, get_requested_port_or_default
from narupa.command.command_service import CommandService
from narupa.protocol.command import add_CommandServicer_to_server

DEFAULT_PORT = 54324


class CommandServer(GrpcServer):
    """
    An implementation of a server for use with :class:`CommandServicer`, primarily for testing.
    """

    def __init__(self, *, address: Optional[str] = None, port: Optional[int] = None):
        if address is None:
            address = DEFAULT_SERVE_ADDRESS
        port = get_requested_port_or_default(port, DEFAULT_PORT)
        super().__init__(address=address, port=port)

    @property
    def service(self) -> CommandService:
        """
        Gets the command service implementation attached to this server.
        :return: The command service.
        """
        return self._service

    def setup_services(self):
        """
        Sets up a new command service and attaches it to the server.
        """
        super().setup_services()
        self._service = CommandService()
        add_CommandServicer_to_server(self._service, self.server)
