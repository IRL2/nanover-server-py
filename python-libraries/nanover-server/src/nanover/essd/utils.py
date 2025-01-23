import ipaddress
import socket
from typing import List, Optional, Iterable, Dict, Union, Literal

import netifaces

Address = str
AddressType = Union[
    Literal["addr"],
    Literal["peer"],
    Literal["mask"],
    Literal["broadcast"],
]

InterfaceAddresses = Dict[AddressType, Address]


def add_broadcast(addr: InterfaceAddresses):
    addr["broadcast"] = str(ipaddress.IPv4Network(f"{addr["addr"]}/{addr["mask"]}", strict=False).broadcast_address)
    return addr


def get_ipv4_addresses(
    interfaces: Optional[Iterable[str]] = None,
) -> List[InterfaceAddresses]:
    """
    Gets all the IPV4 addresses currently available on all the given interfaces.

    :param interfaces: Optional list of interfaces to extract addresses from. If none are provided,
        all interfaces will be used.
    :return: A list of dictionaries containing the IP address and other information for each interface,
        as returned by :func:`netifaces.ifaddresses`.
    """
    if interfaces is None:
        interfaces = netifaces.interfaces()

    ipv4_addrs: List[InterfaceAddresses] = []
    for interface in interfaces:
        addrs = netifaces.ifaddresses(interface)
        try:
            ipv4_addrs += [add_broadcast(addr) for addr in addrs[netifaces.InterfaceType.AF_INET]]
        except KeyError:
            continue

    print(ipv4_addrs)

    return ipv4_addrs


def get_broadcast_addresses(
    interfaces: Optional[Iterable[str]] = None,
) -> List[InterfaceAddresses]:
    """
    Gets all the IPV4 addresses currently available on all the given interfaces that have broadcast addresses.

    :param interfaces: Optional list of interfaces to extract addresses from. If none are provided,
        all interfaces will be used.
    :return: A list of dictionaries containing the IP address and other information for each interface,
        as returned by :func:`netifaces.ifaddresses`.

    In the netifaces API, the address entries are returned as dictionaries in the following format:

    .. code::

        {
          'addr': '172.23.43.33',
          'mask': '255.255.0.0',
          'broadcast': '172.23.255.255'
        }

    """

    ipv4_addrs = get_ipv4_addresses(interfaces)
    return [
        address_entry for address_entry in ipv4_addrs if "broadcast" in address_entry
    ]


def resolve_host_broadcast_address(
    host: str,
    ipv4_addrs: Optional[List[InterfaceAddresses]] = None,
):
    try:
        address = socket.gethostbyname(host)
    except socket.error:
        return None
    if ipv4_addrs is None:
        ipv4_addrs = get_ipv4_addresses()
    return next(
        (
            item
            for item in ipv4_addrs
            if item["addr"] == address and "broadcast" in item
        ),
        None,
    )


def is_in_network(address: str, interface_address_entry: InterfaceAddresses) -> bool:
    """
    An internal mechanism for determining whether a given IP address is part of the same network as a given
    interface network as defined by their IPv4 subnet mask and broadcast address.

    :param address: An IPv4 address.
    :param interface_address_entry: An IPv4 address entry, as produced by :func:`netifaces.ifaddresses`. It must
        contain the `mask` and `broadcast` fields, representing the subnet mask IP and the broadcast IP for the given
        interface
    :return: `True`, if the given address is in the same network as given interface address, `False` otherwise.
    :raises: ValueError: if invalid IP addresses are given for any field.
    :raises: KeyError: if the `mask` and `broadcast` fields are not present in the interface address entry
        argument.
    """
    try:
        ip_address = ipaddress.ip_address(address)
    except ValueError:
        raise ValueError(f"Given address {address} is not a valid IP address.")
    try:
        mask = ipaddress.ip_address(interface_address_entry["mask"])
        broadcast_address = ipaddress.ip_address(interface_address_entry["broadcast"])
        # to network address e.g. 255.255.255.0 & 192.168.1.255 = 192.168.1.0
        network_address = ipaddress.ip_address(int(mask) & int(broadcast_address))
        # The doc and typing stub seem to indicate this is not a valid call of
        # ipaddress.ip_network, but this is well tested so we accept it for the
        # time being.
        # TODO: Fix this line as the types seem to be incorrect.
        ip_network = ipaddress.ip_network((network_address, interface_address_entry["mask"]))  # type: ignore
    except ValueError:
        raise ValueError(
            f"Given address {interface_address_entry} is not a valid IP network address."
        )
    except KeyError:
        raise KeyError(
            f"Given interface address dictionary did not contain either 'broadcast' or 'mask' keys: "
            f"{interface_address_entry}"
        )
    return ip_address in ip_network


def get_broadcastable_ip():
    broadcast_addresses = get_broadcast_addresses()
    if len(broadcast_addresses) == 0:
        raise RuntimeError(
            "No broadcastable IP addresses could be found on the system!"
        )
    return broadcast_addresses[0]["addr"]
