import time

import pytest

from narupa.essd import DiscoveryServer
from narupa.essd.client import DiscoveryClient
from test_essd_server import service
from test_essd_service import properties

from narupa.essd.servicehub import ServiceHub

TEST_SEARCH_TIME = 0.6
TEST_INTERVAL_TIME = 0.001


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


@pytest.mark.timeout(0.4)
def test_client_timeout(client):
    services = client.search_for_services(search_time=0.2)
    assert not services


def test_send_service(client_server, service):
    client, server = client_server
    server.register_service(service)
    services = client.search_for_services(search_time=TEST_SEARCH_TIME, interval=TEST_INTERVAL_TIME)
    assert len(services) == 1
    assert service in services


def test_send_service_different_port(service):
    with DiscoveryServer(broadcast_port=8923) as server:
        with DiscoveryClient(port=8923) as client:
            server.register_service(service)
            services = client.search_for_services(search_time=TEST_SEARCH_TIME, interval=TEST_INTERVAL_TIME)
            assert len(services) == 1
            assert service in services


def test_remove_service(client_server, service):
    client, server = client_server
    server.register_service(service)
    services = client.search_for_services(search_time=TEST_SEARCH_TIME, interval=TEST_INTERVAL_TIME)
    assert len(services) == 1
    assert service in services
    server.unregister_service(service)
    time.sleep(TEST_INTERVAL_TIME)
    services = client.search_for_services(search_time=TEST_SEARCH_TIME, interval=TEST_INTERVAL_TIME)
    assert service not in services


def test_send_service_multiple_clients(client_server, service):
    client, server = client_server
    with DiscoveryClient(port=client.port) as second_client:
        server.register_service(service)
        services = client.search_for_services(search_time=TEST_SEARCH_TIME, interval=TEST_INTERVAL_TIME)
        second_services = second_client.search_for_services(search_time=0.4, interval=TEST_INTERVAL_TIME)
        assert len(services) == 1 == len(second_services)
        assert service in services
        assert service in second_services


def test_send_service_all_interfaces(client_server, service):
    client, server = client_server
    properties = dict(service.properties)
    properties['address'] = '[::]'
    service = ServiceHub(**properties)
    server.register_service(service)
    services = client.search_for_services(search_time=TEST_SEARCH_TIME, interval=TEST_INTERVAL_TIME)
    assert len(services) == 1
    assert ServiceHub(name=service.name, address='127.0.0.1', id=service.id) in services


def run_with_client(service):
    with DiscoveryClient() as client:
        services = client.search_for_services(search_time=TEST_SEARCH_TIME, interval=TEST_INTERVAL_TIME)
        assert len(services) == 1
        assert service in services


def run_with_server(service):
    with DiscoveryServer() as server:
        server.register_service(service)
        for i in range(3):
            run_with_client(service)


@pytest.mark.serial
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
def test_send_utf8(client_server, service, utf_str):
    client, server = client_server
    service.properties['name'] = service.name + utf_str
    server.register_service(service)
    services = client.search_for_services(search_time=TEST_SEARCH_TIME, interval=TEST_INTERVAL_TIME)
    assert len(services) == 1
    assert service in services
