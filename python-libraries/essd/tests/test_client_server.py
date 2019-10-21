import pytest

from essd.client import DiscoveryClient
from test_service import properties
from test_server import server, service

from essd.servicehub import ServiceHub


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


def test_send_service_localhost(client, server, service):
    properties = dict(service.properties)
    properties['address'] = 'localhost'
    service = ServiceHub(**properties)
    server.register_service(service)
    services = client.search_for_services(search_time=0.4, interval=0.001)
    assert len(services) == 1
    assert ServiceHub(name=service.name, address='127.0.0.1') in services


def test_send_service_all_interfaces(client, server, service):
    properties = dict(service.properties)
    properties['address'] = '[::]'
    service = ServiceHub(**properties)
    server.register_service(service)
    services = client.search_for_services(search_time=0.4, interval=0.001)
    assert len(services) == 1
    assert ServiceHub(name=service.name, address='127.0.0.1') in services
