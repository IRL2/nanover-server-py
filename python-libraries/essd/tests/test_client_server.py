import pytest

from essd.client import DiscoveryClient
from test_service import properties
from test_server import server, service


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
