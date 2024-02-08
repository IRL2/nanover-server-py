"""
Tests for application level autoconnecting between client and server.
"""

import pytest
from mock import Mock
from nanover.app import NanoverImdApplication, NanoverImdClient
from nanover.app.app_server import MULTIPLAYER_SERVICE_NAME
from nanover.core import NanoverServer
from nanover.essd import DiscoveryServer, ServiceHub
from nanover.essd.utils import get_broadcastable_ip
from nanover.imd import ImdServer, IMD_SERVICE_NAME
from nanover.trajectory import FrameServer, FRAME_SERVICE_NAME

DISCOVERY_DELAY = 0.05
AUTOCONNECT_SEARCH_TIME = DISCOVERY_DELAY * 1.5

NEVER_USED_HUB_NAME = "pytest adult yoda"


@pytest.fixture
def broadcastable_servers():
    """
    Returns frame server, imd server and multiplayer server configured to be running at
    broadcastable addresses.

    We do this because localhost is not always broadcastable in test environments.
    """
    address = get_broadcastable_ip()
    with FrameServer(address=address, port=0) as frame_server:
        with ImdServer(address=address, port=0) as imd_server:
            with NanoverServer(address=address, port=0) as multiplayer_server:
                yield frame_server, imd_server, multiplayer_server


@pytest.fixture
def discoverable_imd_server():
    """
    Returns a discoverable iMD server discoverable on the free port.
    """
    DISCOVERY_PORT = 39420
    address = get_broadcastable_ip()
    server = NanoverServer(address=address, port=0)
    discovery = DiscoveryServer(broadcast_port=DISCOVERY_PORT, delay=DISCOVERY_DELAY)
    with NanoverImdApplication(server, discovery) as app_server:
        yield app_server


@pytest.fixture
def discoverable_imd_server():
    """
    Returns a discoverable iMD server discoverable on the free port.
    """
    DISCOVERY_PORT = 39421
    address = get_broadcastable_ip()
    server = NanoverServer(address=address, port=0)
    discovery = DiscoveryServer(broadcast_port=DISCOVERY_PORT, delay=DISCOVERY_DELAY)
    with NanoverImdApplication(server, discovery) as app_server:
        yield app_server


@pytest.mark.serial
def test_autoconnect_app_server_default_ports(discoverable_imd_server):
    """
    Tests that an iMD application server running on one port one default port is discoverable and
    that the client connects to it in the expected way.
    """
    mock = Mock(return_value={})

    discoverable_imd_server.server.register_command("test", mock)

    address = get_broadcastable_ip()
    with NanoverImdApplication.basic_server(address=address):
        with NanoverImdClient.autoconnect(
            search_time=AUTOCONNECT_SEARCH_TIME,
            discovery_port=discoverable_imd_server.discovery.port,
        ) as client:
            assert (
                len(client._channels) == 1
            )  # expect the client to connect to each server on the same channel
            # since the command is registered only once on the server, calling it from different 'clients' will
            # actually result in the same method being called 3 times.
            client.run_trajectory_command("test")
            client.run_imd_command("test")
            client.run_multiplayer_command("test")
            assert mock.call_count == 3


def test_autoconnect_app_server(discoverable_imd_server):
    """
    Tests that an iMD application server running on one port is discoverable and
    that the client connects to it in the expected way.
    """
    mock = Mock(return_value={})

    discoverable_imd_server.server.register_command("test", mock)

    with NanoverImdClient.autoconnect(
        search_time=AUTOCONNECT_SEARCH_TIME,
        discovery_port=discoverable_imd_server.discovery.port,
    ) as client:
        assert (
            len(client._channels) == 1
        )  # expect the client to connect to each server on the same channel
        # since the command is registered only once on the server, calling it from different 'clients' will
        # actually result in the same method being called 3 times.
        client.run_trajectory_command("test")
        client.run_imd_command("test")
        client.run_multiplayer_command("test")
        assert mock.call_count == 3


def test_autoconnect_separate_servers(broadcastable_servers):
    """
    Tests that an iMD application running on multiple separate servers on multiple ports is discoverable
    and that the client connects to it in the expected way.
    """
    DISCOVERY_PORT = 39423
    frame_server, imd_server, multiplayer_server = broadcastable_servers

    frame_mock = Mock(return_value={})
    imd_mock = Mock(return_value={})
    multiplayer_mock = Mock(return_value={})

    frame_server.register_command("frame", frame_mock)
    imd_server.register_command("imd", imd_mock)
    multiplayer_server.register_command("multiplayer", multiplayer_mock)

    service_hub = ServiceHub(name="test", address=frame_server.address)
    service_hub.add_service(FRAME_SERVICE_NAME, frame_server.port)
    service_hub.add_service(IMD_SERVICE_NAME, imd_server.port)
    service_hub.add_service(MULTIPLAYER_SERVICE_NAME, multiplayer_server.port)

    with DiscoveryServer(
        broadcast_port=DISCOVERY_PORT, delay=DISCOVERY_DELAY
    ) as discovery_server:
        discovery_server.register_service(service_hub)
        with NanoverImdClient.autoconnect(
            search_time=AUTOCONNECT_SEARCH_TIME, discovery_port=discovery_server.port
        ) as client:
            assert (
                len(client._channels) == 3
            )  # expect the client to connect to each server on a separate channel
            # test servers by running a command on each.
            client.run_trajectory_command("frame")
            client.run_imd_command("imd")
            client.run_multiplayer_command("multiplayer")
            frame_mock.assert_called_once()
            imd_mock.assert_called_once()
            multiplayer_mock.assert_called_once()


@pytest.mark.serial
def test_autoconnect_named_server():
    """
    Test autoconnecting to a named server.
    """
    DISCOVERY_PORT = 39420
    SERVER_NAME = "pytest baby yoda"
    address = get_broadcastable_ip()
    server = NanoverServer(address=address, port=0)
    discovery = DiscoveryServer(broadcast_port=DISCOVERY_PORT, delay=DISCOVERY_DELAY)

    with NanoverImdApplication(server, discovery, name=SERVER_NAME):
        with NanoverImdClient.autoconnect(
            search_time=AUTOCONNECT_SEARCH_TIME,
            discovery_port=DISCOVERY_PORT,
            name=SERVER_NAME,
        ):
            pass


@pytest.mark.serial
def test_autoconnect_no_named_server(discoverable_imd_server):
    """
    Test that autoconnecting to a named server that doesn't exist fails.
    """
    with pytest.raises(ConnectionError), NanoverImdClient.autoconnect(
        name=NEVER_USED_HUB_NAME
    ):
        pass
