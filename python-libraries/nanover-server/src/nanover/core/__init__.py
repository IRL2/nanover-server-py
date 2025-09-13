"""
Module providing core classes and utility functions for NanoVer applications.
"""

from typing import Optional
from .nanover_server import NanoverServer

DEFAULT_SERVE_ADDRESS = "[::]"
DEFAULT_CONNECT_ADDRESS = "localhost"


def get_requested_port_or_default(port: Optional[int], default: int) -> int:
    """
    Returns the port you asked for, or the default one is `port` is `None`.
    """
    if port is None:
        port = default
    return port
