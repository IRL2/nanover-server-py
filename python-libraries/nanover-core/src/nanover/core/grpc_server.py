"""
Module providing a wrapper around the running of GRPC servers.
"""

import logging
from concurrent import futures
from typing import Optional, Tuple

import grpc

DEFAULT_SERVE_ADDRESS = "[::]"
DEFAULT_CONNECT_ADDRESS = "localhost"

# We expect that reserving a large number of threads should not present a
# performance issue. Each concurrent GRPC request requires a worker, and streams
# occupy those workers indefinitely, so several workers must be available for
# each expected client.
DEFAULT_MAX_WORKERS = 128


class GrpcServer:
    """
    A base class for running GRPC servers that handles the starting and closing
    of the underlying server.

    :param address: The IP address at which to run the server.
    :param port: The port on which to run the server.
    """

    def __init__(
        self,
        *,
        address: str,
        port: int,
        max_workers=DEFAULT_MAX_WORKERS,
    ):
        grpc_options = (
            # do not allow hosting two servers on the same port
            ("grpc.so_reuseport", 0),
        )
        executor = futures.ThreadPoolExecutor(max_workers=max_workers)
        self.server = grpc.server(executor, options=grpc_options)
        self.setup_services()
        self._address = address
        try:
            self._port = self.server.add_insecure_port(address=f"{address}:{port}")
        except RuntimeError:
            if port == 0:
                raise IOError("Could not open any port.")
            raise IOError(f"Could not open on port {port}.")
        self._address = address

        self.logger = logging.getLogger(__name__)
        self.logger.info(
            f"Running server {self.__class__.__name__} on port {self.port}"
        )
        self.server.start()

    @property
    def address(self):
        """
        Get the address that this server is or was provided at.
        """
        return self._address

    @property
    def port(self):
        """
        Get the port that the server is or was provided on. This is 0 if a port
        was unable to be chosen.
        """
        return self._port

    @property
    def address_and_port(self) -> Tuple[str, int]:
        """
        Gets the address and port that the server is or was provided on as a tuple.

        :return: The address and port that the server is or was provided on as a tuple.
        """
        return self.address, self.port

    def setup_services(self):
        """
        Inheritors of this class should setup any services they run.
        """
        pass

    def add_service(self, service):
        """
        Add a gRPC service to this server.

        :param service: The gRPC service to add, must have the method to add the gRPC service as the attribute
        add_to_server_method.
        """
        try:
            service.add_to_server_method(service, self.server)
        except AttributeError:
            raise AttributeError(
                "Service implementation did not have the add_to_server_method "
                "as an attribute, cannot automatically add to gRPC server."
            )

    def close(self):
        """
        Stops the server.

        Inheritors of this class should override this method with routines to stop
        services that are running.
        """
        self.server.stop(grace=False)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def get_requested_port_or_default(port: Optional[int], default: int) -> int:
    """
    Returns the port you asked for, or the default one is `port` is `None`.
    """
    if port is None:
        port = default
    return port
