"""
A module containing a Extremely Simple Service Discovery client.
"""

import json
import time
from typing import Optional, Set, Iterable

import select
from socket import socket, AF_INET, SOCK_DGRAM
from _socket import SO_REUSEADDR, SO_REUSEPORT, SOL_SOCKET

from nanover.essd.server import BROADCAST_PORT
from nanover.essd.servicehub import ServiceHub, MAXIMUM_MESSAGE_SIZE

IP_ADDRESS_ANY = "0.0.0.0"


def configure_reusable_socket() -> socket:
    """
    Sets up a socket set up for listening with reuseable address.

    :return: A socket.
    """
    # IPv4 UDP socket
    s = socket(AF_INET, SOCK_DGRAM)
    # Enable reuse
    s.setsockopt(SOL_SOCKET, SO_REUSEADDR, 1)
    s.setsockopt(SOL_SOCKET, SO_REUSEPORT, 1)
    return s


class DiscoveryClient:
    def __init__(self, address: Optional[str] = None, port: Optional[int] = None):
        if address is None:
            address = IP_ADDRESS_ANY
        if port is None:
            port = BROADCAST_PORT
        self.address = address
        self._connect(port)

    @property
    def port(self):
        return self._socket.getsockname()[1]

    def _connect(self, port):
        self._socket = configure_reusable_socket()
        self._socket.bind((self.address, port))

    def _check_for_messages(self, timeout):
        socket_list = [self._socket]
        readable, _, exceptional = select.select(socket_list, [], socket_list, timeout)
        if len(exceptional) > 0:
            raise ConnectionError("Exception on socket while checking for messages.")
        return len(readable) > 0

    def _receive_service(self):
        (message, address) = self._socket.recvfrom(MAXIMUM_MESSAGE_SIZE)
        properties = json.loads(message.decode())
        return ServiceHub(**properties)

    def search_for_services(
        self, search_time: float = 5.0, interval=0.033
    ) -> Iterable[ServiceHub]:
        """
        Searches for and yields services for the given search time.

        :param search_time: Time, in seconds, to search for.
        :param interval: Interval in seconds to wait between checking for new service broadcasts.
        :return: A set of services discovered over the duration.
        """
        services: Set[ServiceHub] = set()
        deadline = time.monotonic() + search_time
        while time.monotonic() < deadline:
            time_before_recv = time.monotonic()
            time_remaining = deadline - time.monotonic()
            if self._check_for_messages(timeout=time_remaining):
                service = self._receive_service()
                if service is not None and service not in services:
                    services.add(service)
                    yield service
            time_spent_receiving = time.monotonic() - time_before_recv
            time_remaining = interval - time_spent_receiving
            if time_remaining > 0:
                time.sleep(time_remaining)

    def close(self):
        self._socket.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
