import ipaddress
import socket
from typing import List, Optional, Iterable, Dict, Any

import netifaces

InterfaceAddresses = Dict[str, str]


def get_ipv4_addresses(interfaces: Optional[Iterable[str]] = None) -> List[InterfaceAddresses]:
    """
    Gets all the IPV4 addresses currently available on all the given interfaces.

    :param interfaces: Optional list of interfaces to extract addresses from. If none are provided,
        all interfaces will be used.
    :return: A list of dictionaries containing the IP address and other information for each interface,
        as returned by :fun:`netifaces.ifaddresses`.
    """
    if interfaces is None:
        interfaces = netifaces.interfaces()
    ipv4_addrs: List[InterfaceAddresses] = []
    for interface in interfaces:
        addrs = netifaces.ifaddresses(interface)
        try:
            ipv4_addrs += addrs[netifaces.AF_INET]
        except KeyError:
            continue
    return ipv4_addrs


def get_broadcast_addresses(interfaces: Optional[Iterable[str]] = None) -> List[Dict[str, str]]:
    """
    Gets all the IPV4 addresses currently available on all the given interfaces that have broadcast addresses.

    :param interfaces: Optional list of interfaces to extract addresses from. If none are provided,
        all interfaces will be used.
    :return: A list of dictionaries containing the IP address and other information for each interface,
    as returned by :fun:`netifaces.ifaddresses`.

    In the netifaces API, the address entries are returned as dictionaries in the following format:

    .. code::
        {
          'addr': '172.23.43.33',
          'netmask': '255.255.0.0',
          'broadcast': '172.23.255.255'
        }

    """

    ipv4_addrs = get_ipv4_addresses(interfaces)
    return [address_entry for address_entry in ipv4_addrs if 'broadcast' in address_entry]


def resolve_host_broadcast_address(
        host: str,
        ipv4_addrs: List[InterfaceAddresses] = None,
):
    try:
        address = socket.gethostbyname(host)
    except socket.error:
        return None
    if ipv4_addrs is None:
        ipv4_addrs = get_ipv4_addresses()
    return next((item for item in ipv4_addrs if item["addr"] == address and 'broadcast' in item), None)


def is_in_network(address: str, interface_address_entry: InterfaceAddresses) -> bool:
    """
    An internal mechanism for determining whether a given IP address is part of the same network as a given
    interface network as defined by their IPv4 subnet mask and broadcast address.

    :param address: An IPv4 address.
    :param interface_address_entry: An IPv4 address entry, as produced by :fun:`netifaces.ifaddresses`. It must
        contain the `netmask` and `broadcast` fields, representing the subnet mask IP and the broadcast IP for the given
        interface
    :return: `True`, if the given address is in the same network as given interface address, `False` otherwise.
    :raises: ValueError: if invalid IP addresses are given for any field.
    :raises: KeyError: if the `netmask` and `broadcast` fields are not present in the interface address entry
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


def get_broadcastable_ip():
    broadcast_addresses = get_broadcast_addresses()
    if len(broadcast_addresses) == 0:
        raise RuntimeError("No broadcastable IP addresses could be found on the system!")
    return broadcast_addresses[0]['addr']