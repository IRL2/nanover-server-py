# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing service discovery server.

"""
import ipaddress
import logging
import threading
import time
from collections import namedtuple
from socket import socket, AF_INET, SOCK_DGRAM, SOL_SOCKET, SO_BROADCAST
from typing import Optional, List

import netifaces

from .servicehub import ServiceHub

BROADCAST_PORT = 54545
IP_ADDRESS_BROADCAST = '255.255.255.255'

ServiceHubEntry = namedtuple('ServiceHubEntry', ['service', 'broadcast_addresses'])


def _connect_socket() -> socket:
    # IPv4 UDP socket
    s = socket(AF_INET, SOCK_DGRAM)
    # Enable broadcasting
    s.setsockopt(SOL_SOCKET, SO_BROADCAST, 1)
    return s


def get_ipv4_addresses(interfaces: List[str] = None):
    """
    Gets all the IPV4 addresses currently available on all the given interfaces.
    :param interfaces: Optional list of interfaces to extract addresses from. If none are provided,
    all interfaces will be used.
    :return: A list of dictionaries containing the IP address and other information for each interface,
    as returned by :fun:`netifaces.ifaddresses`.
    """
    if interfaces is None:
        interfaces = netifaces.interfaces()
    ipv4_addrs = []
    for interface in interfaces:
        addrs = netifaces.ifaddresses(interface)
        try:
            ipv4_addrs += addrs[netifaces.AF_INET]
        except KeyError:
            continue
    return ipv4_addrs


def get_broadcast_addresses(interfaces: List[str] = None):
    """
    Gets all the IPV4 addresses currently available on all the given interfaces that have broadcast addresses.
    :param interfaces: Optional list of interfaces to extract addresses from. If none are provided,
    all interfaces will be used.
    :return: A list of dictionaries containing the IP address and other information for each interface,
    as returned by :fun:`netifaces.ifaddresses`.

    In the netifaces API, the address entries are returned as dictionaries in the following format:

    .. code
        {
          'addr': '172.23.43.33',
          'netmask': '255.255.0.0',
          'broadcast': '172.23.255.255'
        }

    """

    ipv4_addrs = get_ipv4_addresses(interfaces)
    return [address_entry for address_entry in ipv4_addrs if 'broadcast' in address_entry]


def is_in_network(address, interface_address_entry):
    """
    An internal mechanism for determining whether a given IP address is part of the same network as a given
    interface network as defined by their IPv4 subnet mask and broadcast address.

    :param address: An IPv4 address.
    :param interface_address_entry: An IPv4 address entry, as produced by :fun:`netifaces.ifaddresses`. It must
    contain the `netmask` and `broadcast` fields, representing the subnet mask IP and the broadcast IP for the given
    interface.
    :return: `True`, if the given address is in the same network as given interface address, `False` otherwise.
    :raises: ValueError, if invalid IP addresses are given for any field.
    :raises: KeyError, if the `netmask` and `broadcast` fields are not present in the interface address entry
    argument.
    """
    try:
        ip_address = ipaddress.ip_address(address)
    except ValueError:
        raise ValueError(f'Given address {address} is not a valid IP address.')
    try:
        netmask = ipaddress.ip_address(interface_address_entry['netmask'])
        broadcast_address = ipaddress.ip_address(interface_address_entry['broadcast'])
        # to network address e.g. 255.255.255.0 & 192.168.1.255 = 192.168.1.0
        network_address = ipaddress.ip_address(int(netmask) & int(broadcast_address))
        ip_network = ipaddress.ip_network((network_address, interface_address_entry['netmask']))
    except ValueError:
        raise ValueError(f'Given address {interface_address_entry} is not a valid IP network address.')
    except KeyError:
        raise KeyError(f'Given interface address dictionary did not contain either \'broadcast\' or \'netmask\' keys: '
                       f'{interface_address_entry}')
    return ip_address in ip_network


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
        self.logger.log(level, f"ESSD: Able to broadcast on the following IPV4 addresses:")
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
            raise ValueError(f"No valid broadcast address found for service {service}, check network configuration.")
        with self._lock:
            self.services[service] = broadcast_addresses

    @property
    def broadcasting(self):
        return self._broadcast_thread is not None

    def start(self):
        if self._broadcast_thread is not None:
            raise RuntimeError("Discovery service already running!")
        self._socket = _connect_socket()
        self._broadcast_thread = threading.Thread(target=self._broadcast, daemon=True)
        self._broadcast_thread.start()

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
            for service, addresses in self.services.items():
                self._broadcast_service(service, addresses)

    def _broadcast_service(self, service: ServiceHub, addresses):
        address = service.address
        for broadcast_address in addresses:
            if address == "[::]" or address == "localhost":
                message = service.to_message(override_address=broadcast_address['addr'])
            else:
                message = service.to_message()
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
            # manually construct an address/broadcast address pair for the localhost shortcut.
            return [{'addr': '127.0.0.1', 'broadcast': '127.255.255.255'}]
        return [broadcast_address for broadcast_address in self.broadcast_addresses
                if is_in_network(address, broadcast_address)]
