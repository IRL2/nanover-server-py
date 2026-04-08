import ipaddress
import socket

import psutil
from psutil._common import snicaddr


def snicaddr_with_computed_broadcast_address(addr: snicaddr):
    if addr.broadcast is None and addr.ptp is None:
        broadcast = str(
            ipaddress.IPv4Network(
                f"{addr.address}/{addr.netmask}", strict=False
            ).broadcast_address
        )
        return addr._replace(broadcast=broadcast)
    else:
        return addr


def get_ipv4_addresses() -> list[snicaddr]:
    """
    Gets all the IPV4 addresses currently available on all interfaces.

    :return: A list of dictionaries containing the IP address and other information for each interface,
        as returned by :func:`netifaces.ifaddresses`.
    """
    active_ifs = {name for name, stats in psutil.net_if_stats().items() if stats.isup}
    valid_ifs = {
        name: addrs
        for name, addrs in psutil.net_if_addrs().items()
        if name in active_ifs
    }

    ipv4_addrs = [
        snicaddr_with_computed_broadcast_address(addr)
        for name, addrs in valid_ifs.items()
        for addr in addrs
        if addr.family == socket.AddressFamily.AF_INET
    ]

    return ipv4_addrs


def get_broadcast_addresses() -> list[snicaddr]:
    """
    Gets all the IPV4 addresses currently available on all interfaces that have broadcast addresses.

    :return: A list of dictionaries containing the IP address and other information for each interface,
        as returned by :func:`netifaces.ifaddresses`.

    In the netifaces API, the address entries are returned as dictionaries in the following format:

    .. code::

        {
          'addr': '172.23.43.33',
          'netmask': '255.255.0.0',
          'broadcast': '172.23.255.255'
        }

    """

    ipv4_addrs = get_ipv4_addresses()
    print(ipv4_addrs)
    return [
        address_entry
        for address_entry in ipv4_addrs
        if address_entry.broadcast is not None
    ]


def resolve_host_broadcast_address(
    host: str,
    ipv4_addrs: list[snicaddr] | None = None,
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
            if item.address == address and item.broadcast is not None
        ),
        None,
    )


def is_in_network(address: str, interface_address_entry: snicaddr) -> bool:
    """
    An internal mechanism for determining whether a given IP address is part of the same network as a given
    interface network as defined by their IPv4 subnet mask and broadcast address.

    :param address: An IPv4 address.
    :param interface_address_entry: An IPv4 address entry, as produced by :func:`netifaces.ifaddresses`. It must
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
        raise ValueError(f"Given address {address} is not a valid IP address.")
    try:
        netmask = ipaddress.ip_address(interface_address_entry.netmask)
        broadcast_address = ipaddress.ip_address(interface_address_entry.broadcast)
        # to network address e.g. 255.255.255.0 & 192.168.1.255 = 192.168.1.0
        network_address = ipaddress.ip_address(int(netmask) & int(broadcast_address))
        # The doc and typing stub seem to indicate this is not a valid call of
        # ipaddress.ip_network, but this is well tested so we accept it for the
        # time being.
        # TODO: Fix this line as the types seem to be incorrect.
        ip_network = ipaddress.ip_network(
            (network_address, interface_address_entry.netmask)  # type: ignore
        )
    except ValueError:
        raise ValueError(
            f"Given address {interface_address_entry} is not a valid IP network address."
        )
    except KeyError:
        raise KeyError(
            f"Given interface address dictionary did not contain either 'broadcast' or 'netmask' keys: "
            f"{interface_address_entry}"
        )
    return ip_address in ip_network


def get_broadcastable_test_ip():
    broadcast_addresses = get_broadcast_addresses()
    if len(broadcast_addresses) == 0:
        raise RuntimeError(
            "No broadcastable IP addresses could be found on the system!"
        )

    return broadcast_addresses[0].address
