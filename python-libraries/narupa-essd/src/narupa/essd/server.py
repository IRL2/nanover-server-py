# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing a server for Extremely Simple Service Discovery (ESSD).

Narupa servers can use this class to broadcast themselves on a local area network, using UDP.

Example
=======

>>> discovery_server = DiscoveryServer(broadcast_port=54545)
>>> # assume one has created a multiplayer server.
>>> hub = ServiceHub(name="Example Narupa Service Hub", address="[::]")
>>> hub.add_service("multiplayer", 54323)
>>> discovery_server.close()

"""
import logging
import threading
import time
from socket import socket, AF_INET, SOCK_DGRAM, SOL_SOCKET, SO_BROADCAST, SO_REUSEADDR
from typing import Optional

from narupa.essd.utils import get_broadcast_addresses, is_in_network, resolve_host_broadcast_address
from narupa.essd.servicehub import ServiceHub

BROADCAST_PORT = 54545


def configure_reusable_socket() -> socket:
    """
    Sets up a socket set up for broadcasting with reuseable address.

    :return: A socket.
    """
    # IPv4 UDP socket
    s = socket(AF_INET, SOCK_DGRAM)
    # Enable broadcasting
    s.setsockopt(SOL_SOCKET, SO_BROADCAST, 1)
    s.setsockopt(SOL_SOCKET, SO_REUSEADDR, 1)
    return s


class DiscoveryServer:

    def __init__(self, broadcast_port: Optional[int] = None, delay=0.5):
        if broadcast_port is None:
            broadcast_port = BROADCAST_PORT
        self.logger = logging.getLogger(__name__)
        self.port = broadcast_port
        self.logger.info(f"Extremely Simple Discovery Server (ESSD) is running, "
                         f"will broadcast services to port {self.port}")
        self.broadcast_addresses = get_broadcast_addresses()
        self.log_addresses(level=logging.INFO)
        self.delay = delay
        self.services = dict()
        self._lock = threading.RLock()
        self._cancel = False
        self._broadcast_thread = None
        self._socket = None
        self.start()

    def log_addresses(self, level=logging.DEBUG):
        self.logger.log(level, "ESSD: Able to broadcast on the following IPV4 addresses:")
        for address in self.broadcast_addresses:
            self.logger.log(level, f"ESSD:   - {address}")

    def register_service(self, service: ServiceHub):
        """
        Register a service for discovery.

        :param service: Service to register.
        """
        if service in self.services:
            raise KeyError(f"A service with the same ID has already been registered: {service}")
        broadcast_addresses = self.get_broadcast_addresses_for_service(service)
        if len(broadcast_addresses) == 0:
            msg = f"No valid broadcast address found for service {service}, check network configuration. "
            if len(self.broadcast_addresses) > 0:
                msg += f"The following broadcast addresses were found on the system: {self.broadcast_addresses}"
            raise ValueError(msg)
        with self._lock:
            self.services[service] = broadcast_addresses

    def unregister_service(self, service: ServiceHub):
        """
        Removes a service from discovery.

        :param service: The service to remove.
        :raises KeyError Raised if the service has never been registered with this discovery server.
        """
        if service not in self.services:
            raise KeyError(f'No service with this ID has been registered {service}')
        with self._lock:
            del self.services[service]

    @property
    def is_broadcasting(self):
        return self._broadcast_thread is not None

    def start(self):
        if self._broadcast_thread is not None:
            raise RuntimeError("Discovery service already running!")
        self._socket = configure_reusable_socket()
        self._broadcast_thread = threading.Thread(target=self._broadcast, daemon=True)
        self._broadcast_thread.start()

    def close(self):
        if self.is_broadcasting:
            self._cancel = True
            self._broadcast_thread.join()
            self._broadcast_thread = None
            self._cancel = False
            self._socket.close()

    def _broadcast(self):
        while not self._cancel:
            self._broadcast_services()
            time.sleep(self.delay)

    def _broadcast_services(self):
        with self._lock:
            for service, addresses in self.services.items():
                self._broadcast_one_service(service, addresses)

    def _broadcast_one_service(self, service: ServiceHub, addresses):
        address = service.address
        for broadcast_address in addresses:
            if address == "[::]" or address == "localhost":
                message = service.to_message(override_address=broadcast_address['addr'])
            else:
                message = service.to_message()
            self.logger.debug(f'Sending service {service} to {broadcast_address["broadcast"]}:{self.port}')
            self._socket.sendto(message.encode(), (broadcast_address['broadcast'], self.port))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def get_broadcast_addresses_for_service(self, service):
        address = service.address
        if address == "[::]":
            return self.broadcast_addresses
        if address == "localhost":
            localhost_address = resolve_host_broadcast_address(address)
            if localhost_address is None:
                raise ValueError("Cannot broadcast on localhost on this system!")
            return [localhost_address]
        return [broadcast_address for broadcast_address in self.broadcast_addresses
                if is_in_network(address, broadcast_address)]
