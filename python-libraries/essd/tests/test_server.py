import pytest

from essd.server import DiscoveryServer
from essd.servicehub import ServiceHub
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
