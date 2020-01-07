"""
Module providing an out-of-the-box Narupa application server,
with an underyling gRPC server, discovery, multiplayer and commands.
"""
from narupa.core import NarupaServer
from narupa.essd import DiscoveryServer
from narupa.multiplayer.multiplayer_service import MultiplayerService


class NarupaApplicationServer:

    def __init__(self, address=None, port=0):
        self._server = NarupaServer(address=address, port=port)
        self._discovery = DiscoveryServer()
        self._services = []

    def setup_services(self):
        self._multiplayer = MultiplayerService()


    def close(self):
        self._discovery.close()
        for service in self._services:
            service.close()

