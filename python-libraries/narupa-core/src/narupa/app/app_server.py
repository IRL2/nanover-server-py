"""
Module providing an out-of-the-box Narupa application server,
with an underyling gRPC server, discovery, multiplayer and commands.
"""
from narupa.core import NarupaServer
from narupa.essd import DiscoveryServer
from narupa.multiplayer.multiplayer_service import MultiplayerService
from narupa.protocol.multiplayer import add_MultiplayerServicer_to_server


class NarupaApplicationServer:

    def __init__(self, address=None, port=0):
        self._server = NarupaServer(address=address, port=port)
        self._discovery = DiscoveryServer()
        self._services = set()
        self._update_discovery_services()

    def _setup_multiplayer(self):
        self._multiplayer = MultiplayerService()
        self._server.add_service(self._multiplayer, add_MultiplayerServicer_to_server)
        self._services.add(self._multiplayer)
    
    
    def close(self):
        self._discovery.close()
        for service in self._services:
            service.close()
        self._server.close()

    def _update_discovery_services(self):
        pass
