import pytest

from essd.server import DiscoveryServer, get_ipv4_addresses, get_broadcast_addresses, is_in_network
from essd.servicehub import ServiceHub
import netifaces
from test_service import properties


@pytest.fixture
def server():
    server = DiscoveryServer()
    yield server
    server.close()


@pytest.fixture
def service(properties):
    return ServiceHub(**properties)


def test_server(server, service):
    server.register_service(service)


def test_get_ipv4_addresses():
    ipv4_addresses = get_ipv4_addresses()
    assert len(ipv4_addresses) > 0


def test_get_ipv4_addresses_interface():
    interfaces = netifaces.interfaces()
    if interfaces is None or len(interfaces) == 0:
        return
    interface = [interfaces[0]]
    ipv4_addresses = get_ipv4_addresses(interface)
    expected_addresses = netifaces.ifaddresses(interface[0])[netifaces.AF_INET]
    assert ipv4_addresses == expected_addresses


def test_get_broadcast_addresses():
    interfaces = netifaces.interfaces()
    if interfaces is None or len(interfaces) == 0:
        return
    interface = [interfaces[0]]
    broadcast_addresses = get_broadcast_addresses(interface)
    expected_addresses = netifaces.ifaddresses(interface[0])[netifaces.AF_INET]
    assert broadcast_addresses == expected_addresses


@pytest.mark.parametrize('address, netmask, broadcast_address, expected_result',
                         [('192.168.1.2', '255.255.0.0', '192.168.255.255', True),
                          ('192.5.1.2', '255.255.0.0', '192.168.255.255', False),
                          ('192.168.2.3', '255.255.255.0', '10.0.3.255', False),
                          ('192.168.1.2', '255.255.255.0', '192.168.255.255', False),
                          ('127.0.0.1', '255.0.0.0', '127.255.255.255', True),
                          ('127.2.3.4', '255.0.0.0', '127.255.255.255', True)])
def test_is_in_network(address, netmask, broadcast_address, expected_result):
    network_interface_addresses = {'netmask': netmask, 'broadcast': broadcast_address}
    assert expected_result == is_in_network(address, network_interface_addresses)
