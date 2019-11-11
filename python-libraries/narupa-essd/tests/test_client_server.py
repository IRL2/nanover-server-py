import pytest

from narupa.essd import DiscoveryServer
from narupa.essd.client import DiscoveryClient
from test_service import properties
from test_essd_server import server, service

from narupa.essd.servicehub import ServiceHub


@pytest.fixture
def client():
    client = DiscoveryClient()
    yield client
    client.close()


def test_send_service(client, server, service):
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


def test_send_service_multiple_clients(client, server, service):
    with DiscoveryClient() as second_client:
        server.register_service(service)
        services = client.search_for_services(search_time=0.4, interval=0.001)
        second_services = second_client.search_for_services(search_time=0.4, interval=0.001)
        assert len(services) == 1 == len(second_services)
        assert service in services
        assert service in second_services


def test_send_service_all_interfaces(client, server, service):
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
    for i in range(2):
        run_with_server(service)

@pytest.mark.parametrize('utf_str',
                         ['한국어',
                          '😀',
                         ])
def test_send_utf8(client, server, service, utf_str):
    service.properties['name'] = service.name + utf_str
    server.register_service(service)
    services = client.search_for_services(search_time=0.4, interval=0.001)
    assert len(services) == 1
    assert service in services

