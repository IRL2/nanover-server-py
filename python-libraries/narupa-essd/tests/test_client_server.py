import pytest

from narupa.essd import DiscoveryServer
from narupa.essd.client import DiscoveryClient
from test_essd_server import service
from test_service import properties

from narupa.essd.servicehub import ServiceHub


@pytest.fixture
def client():
    client = DiscoveryClient(port=0)
    yield client
    client.close()


@pytest.fixture
def client_server(client):
    server = DiscoveryServer(broadcast_port=client.port)
    yield client, server
    server.close()


def test_send_service(client_server, service):
    client, server = client_server
    server.register_service(service)
    services = client.search_for_services(search_time=0.4, interval=0.001)
    assert len(services) == 1
    assert service in services


def test_send_service_different_port(service):
    with DiscoveryServer(broadcast_port=8923) as server:
        with DiscoveryClient(port=8923) as client:
            server.register_service(service)
            services = client.search_for_services(search_time=0.4, interval=0.001)
            assert len(services) == 1
            assert service in services


def test_send_service_multiple_clients(client_server, service):
    client, server = client_server
    with DiscoveryClient(port=client.port) as second_client:
        server.register_service(service)
        services = client.search_for_services(search_time=0.4, interval=0.001)
        second_services = second_client.search_for_services(search_time=0.4, interval=0.001)
        assert len(services) == 1 == len(second_services)
        assert service in services
        assert service in second_services


def test_send_service_all_interfaces(client_server, service):
    client, server = client_server
    properties = dict(service.properties)
    properties['address'] = '[::]'
    service = ServiceHub(**properties)
    server.register_service(service)
    services = client.search_for_services(search_time=0.4, interval=0.001)
    assert len(services) == 1
    assert ServiceHub(name=service.name, address='127.0.0.1', id=service.id) in services


def run_with_client(service):
    with DiscoveryClient() as client:
        services = client.search_for_services(search_time=0.4, interval=0.001)
        assert len(services) == 1
        assert service in services


def run_with_server(service):
    with DiscoveryServer() as server:
        server.register_service(service)
        for i in range(3):
            run_with_client(service)


def test_context_managers(service):
    """
    tests that running the server and client with context managers cleans up correctly.
    If discovery servers do not clean up cleanly, future clients will find additional servers.
    """
    for port in range(2):
        run_with_server(service)


@pytest.mark.parametrize('utf_str',
                         ['한국어',
                          '😀',
                          ])
def test_send_utf8(client_server, service, utf_str):
    client, server = client_server
    service.properties['name'] = service.name + utf_str
    server.register_service(service)
    services = client.search_for_services(search_time=0.4, interval=0.001)
    assert len(services) == 1
    assert service in services
