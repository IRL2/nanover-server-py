"""
Module containing a basic interactive molecular dynamics client that receives frames
and can publish interactions.
"""

from typing import Optional

from nanover.essd import DiscoveryClient
from nanover.trajectory import FRAME_SERVICE_NAME

# Default to a low framerate to avoid build up in the frame stream
DEFAULT_SUBSCRIPTION_INTERVAL = 1 / 30

# ID of the root selection
SELECTION_ROOT_ID = "selection.root"
# Name of the root selection
SELECTION_ROOT_NAME = "Root Selection"


def _search_for_first_server_with_name(
    server_name: str,
    search_time: float = 2.0,
    discovery_address: Optional[str] = None,
    discovery_port: Optional[int] = None,
):
    with DiscoveryClient(discovery_address, discovery_port) as discovery_client:
        for hub in discovery_client.search_for_services(search_time):
            if hub.name == server_name:
                return hub
    return None


def _search_for_first_available_frame_service(
    search_time: float = 2.0,
    discovery_address: Optional[str] = None,
    discovery_port: Optional[int] = None,
):
    with DiscoveryClient(discovery_address, discovery_port) as discovery_client:
        for hub in discovery_client.search_for_services(search_time):
            if FRAME_SERVICE_NAME in hub.services:
                return hub
    return None
