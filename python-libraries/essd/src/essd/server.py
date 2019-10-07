# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing service discovery server.

"""
import threading
import time
from socket import socket, AF_INET, SOCK_DGRAM, SOL_SOCKET, SO_BROADCAST
from typing import Optional

import essd
from .servicehub import ServiceHub

BROADCAST_PORT = 54545
IP_ADDRESS_BROADCAST = '255.255.255.255'


def get_essd_version() -> str:
    return essd.__version__


def _connect_socket() -> socket:
    # IPv4 UDP socket
    s = socket(AF_INET, SOCK_DGRAM)
    # Enable broadcasting
    s.setsockopt(SOL_SOCKET, SO_BROADCAST, 1)
    return s


class DiscoveryServer:

    def __init__(self, broadcast_address: Optional[str] = None, broadcast_port: Optional[int] = None, delay=0.5):
        if broadcast_address is None:
            broadcast_address = IP_ADDRESS_BROADCAST
        if broadcast_port is None:
            broadcast_port = BROADCAST_PORT

        self.address = broadcast_address
        self.port = broadcast_port
        self.delay = delay
        self.services = set()
        self._lock = threading.RLock()
        self._cancel = False
        self._broadcast_thread = None
        self._socket = None
        self.start()

    def register_service(self, service: ServiceHub):
        """
        Register a service for discovery.
        :param service: Service to register.
        """
        if service in self.services:
            raise KeyError(f"A service with the same name and IP address has already been registered: {service}")
        with self._lock:
            self.services.add(service)

    @property
    def broadcasting(self):
        return self._broadcast_thread is not None

    def start(self):
        if self._broadcast_thread is not None:
            raise RuntimeError("Discovery service already running!")
        self._socket = _connect_socket()
        try:
            self._broadcast_thread = threading.Thread(target=self._broadcast, daemon=True)
            self._broadcast_thread.start()
        except:
            # if anything goes wrong, be sure to close the socket.
            self._socket.close()

    def close(self):
        if self.broadcasting:
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
            for service in self.services:
                self._broadcast_service(service)

    def _broadcast_service(self, service: ServiceHub):
        message = service.message
        self._socket.sendto(message.encode(), (IP_ADDRESS_BROADCAST, BROADCAST_PORT))

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
