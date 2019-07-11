# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing a wrapper around the running of GRPC servers.
"""
from concurrent import futures
from typing import Optional

import grpc

DEFAULT_SERVE_ADDRESS = '[::]'
DEFAULT_CONNECT_ADDRESS = 'localhost'


class GrpcServer:
    """
    A base class for running GRPC servers that handles the starting and closing of the underlying server.

    :param address: The IP address at which to run the server.
    :param port: The port on which to run the server.
    """

    def __init__(self, *, address: str, port: int):
        self.server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))

        self.setup_services()

        self._port = self.server.add_insecure_port(address="{0}:{1}".format(address, port))
        self.server.start()

    def setup_services(self):
        pass

    def close(self):
        self.server.stop(grace=False)


def get_requested_port_or_default(port: Optional[int], default: int) -> int:
    """
    Returns the port you asked for, or the default one is `port` is `None`.
    """
    if port is None:
        port = default
    return port