"""
Tests for application level autoconnecting between client and server.
"""

import pytest
from mock import Mock
from nanover.app import NanoverImdApplication
from nanover.core import AppServer
from nanover.essd import DiscoveryServer
from nanover.essd.server import BROADCAST_PORT
from nanover.essd.utils import get_broadcastable_ip
from nanover.websocket import NanoverImdClient

DISCOVERY_DELAY = 0.05
AUTOCONNECT_SEARCH_TIME = 0.5

NEVER_USED_HUB_NAME = "pytest adult yoda"


@pytest.fixture
def discoverable_imd_server():
    """
    Returns a discoverable iMD server on a free port, discoverable on a non-default ESSD port.
    """
    # Use unique non-default port for discovery. This avoids interference
    # with other tests and other servers on the network.
    DISCOVERY_PORT = BROADCAST_PORT + 1
    address = get_broadcastable_ip()
    discovery = DiscoveryServer(broadcast_port=DISCOVERY_PORT, delay=DISCOVERY_DELAY)
    with NanoverImdApplication(address=address, discovery=discovery) as app_server:
        app_server.serve_websocket()
        yield app_server


@pytest.mark.serial
def test_autoconnect_app_server_default_ports():
    """
    Tests that an iMD application server running on the default ports is discoverable and
    that the client connects to it in the expected way.
    """
    mock = Mock(return_value={})

    address = get_broadcastable_ip()
    discovery = DiscoveryServer(delay=DISCOVERY_DELAY)

    with NanoverImdApplication(discovery=discovery, address=address) as app_server:
        app_server.serve_websocket()
        app_server.register_command("test", mock)
        with NanoverImdClient.from_discovery() as client:
            client.run_command_blocking("test")
            assert mock.call_count == 1


def test_autoconnect_app_server(discoverable_imd_server: AppServer):
    """
    Tests that an iMD application server running on one port is discoverable and
    that the client connects to it in the expected way.
    """
    mock = Mock(return_value={})

    discoverable_imd_server.register_command("test", mock)

    with NanoverImdClient.from_discovery(
        discovery_port=discoverable_imd_server.discovery.port,
    ) as client:
        client.run_command_blocking("test")
        assert mock.call_count == 1


@pytest.mark.serial
def test_autoconnect_named_server():
    """
    Test autoconnecting to a named server.
    """
    SERVER_NAME = "pytest baby yoda"
    address = get_broadcastable_ip()
    discovery = DiscoveryServer(delay=DISCOVERY_DELAY)

    with NanoverImdApplication(
        address=address, discovery=discovery, name=SERVER_NAME
    ) as app_server:
        app_server.serve_websocket()
        with NanoverImdClient.from_discovery(
            server_name=SERVER_NAME,
        ):
            pass


@pytest.mark.serial
def test_autoconnect_no_named_server(discoverable_imd_server):
    """
    Test that autoconnecting to a named server that doesn't exist fails.
    """
    with pytest.raises(ConnectionError), NanoverImdClient.from_discovery(
        server_name=NEVER_USED_HUB_NAME
    ):
        pass
