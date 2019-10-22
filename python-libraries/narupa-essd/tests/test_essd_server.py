import pytest

from narupa.essd.server import DiscoveryServer
from narupa.essd.utils import get_ipv4_addresses, get_broadcast_addresses, is_in_network
from narupa.essd.servicehub import ServiceHub
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


def test_server_duplicate_service(server, service):
    server.register_service(service)
    service_2 = ServiceHub(**service.properties)
    with pytest.raises(KeyError):
        server.register_service(service_2)


def test_server_discovery_already_running(server):
    with pytest.raises(RuntimeError):
        server.start()


def test_server_discovery_restart(server, service):
    server.register_service(service)
    server.close()
    server.start()
    assert service in server.services


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
    broadcast_addresses = get_broadcast_addresses()
    print(broadcast_addresses)
    assert len(broadcast_addresses) > 0



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


@pytest.mark.parametrize('address, netmask, broadcast_address',
                         [('192.168.1.x', '255.255.0.0', '192.168.255.255'),
                          ('192.168.1.2', '255.255.x', '192.168.255.255'),
                          ('192.168.1.2', '255.255.255.0', '192.168.xx.255'),
                          ('192.168.1.2', '255.255.255.0', '192.168.xx.255')])
def test_is_in_network_invalid_addresses(address, netmask, broadcast_address):
    network_interface_addresses = {'netmask': netmask, 'broadcast': broadcast_address}
    with pytest.raises(ValueError):
        _ = is_in_network(address, network_interface_addresses)


@pytest.mark.parametrize('entry',
                         [({'broadcast': '192.168.255.255'}),
                          ({'netmask': '255.255.255.0'}),
                          ({})])
def test_is_in_network_missing_fields(entry):
    with pytest.raises(KeyError):
        _ = is_in_network('192.168.0.1', entry)
