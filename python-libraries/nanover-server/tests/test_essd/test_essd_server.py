import pytest
from psutil._common import snicaddr

from nanover.essd.server import DiscoveryServer
from nanover.essd.utils import (
    get_ipv4_addresses,
    get_broadcast_addresses,
    is_in_network,
    resolve_host_broadcast_address,
)
from nanover.essd.servicehub import ServiceHub
from test_essd_service import (
    properties,
    get_broadcastable_test_ip,
    properties_unique_id,
)


@pytest.fixture
def server():
    server = DiscoveryServer()
    yield server
    server.close()


@pytest.fixture
def service(properties_unique_id):
    return ServiceHub(**properties_unique_id)


def test_server(server, service):
    server.register_service(service)


def test_server_duplicate_service(server, service):
    server.register_service(service)
    service_2 = ServiceHub(**service.properties)
    with pytest.raises(KeyError):
        server.register_service(service_2)


def test_remove_service(server, service):
    server.register_service(service)
    assert service in server.services
    server.unregister_service(service)
    assert service not in server.services


def test_remove_unknown_service(server, service):
    with pytest.raises(KeyError):
        server.unregister_service(service)


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


def test_get_broadcast_addresses():
    broadcast_addresses = get_broadcast_addresses()
    print(broadcast_addresses)
    assert len(broadcast_addresses) > 0


@pytest.mark.parametrize(
    "address, netmask, broadcast_address, expected_result",
    [
        ("192.168.1.2", "255.255.0.0", "192.168.255.255", True),
        ("192.5.1.2", "255.255.0.0", "192.168.255.255", False),
        ("192.168.2.3", "255.255.255.0", "10.0.3.255", False),
        ("192.168.1.2", "255.255.255.0", "192.168.255.255", False),
        ("127.0.0.1", "255.0.0.0", "127.255.255.255", True),
        ("127.2.3.4", "255.0.0.0", "127.255.255.255", True),
    ],
)
def test_is_in_network(address, netmask, broadcast_address, expected_result):
    network_interface_addresses = snicaddr(
        netmask=netmask,
        broadcast=broadcast_address,
        ptp=None,
        address=None,
        family=None,
    )
    assert expected_result == is_in_network(address, network_interface_addresses)


@pytest.mark.parametrize(
    "address, netmask, broadcast_address",
    [
        ("192.168.1.x", "255.255.0.0", "192.168.255.255"),
        ("192.168.1.2", "255.255.x", "192.168.255.255"),
        ("192.168.1.2", "255.255.255.0", "192.168.xx.255"),
        ("192.168.1.2", "255.255.255.0", "192.168.xx.255"),
    ],
)
def test_is_in_network_invalid_addresses(address, netmask, broadcast_address):
    network_interface_addresses = snicaddr(
        netmask=netmask,
        broadcast=broadcast_address,
        ptp=None,
        address=None,
        family=None,
    )
    with pytest.raises(ValueError):
        _ = is_in_network(address, network_interface_addresses)


def test_resolve_address():
    """
    Tests that we can resolve a broadcast address, given a valid address on the network.
    The resolve address function is primarily used with 'localhost', but that
    does not exist on the CI, so we test what we can.
    """
    ip = get_broadcastable_test_ip()
    addr = resolve_host_broadcast_address(ip)
    assert addr is not None
    assert addr.address == ip
    assert addr.broadcast is not None


def test_resolve_invalid_address():
    addr = resolve_host_broadcast_address("blah")
    assert addr is None
