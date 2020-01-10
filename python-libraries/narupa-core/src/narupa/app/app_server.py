"""
Module providing an out-of-the-box Narupa application server,
with an underyling gRPC server, discovery, multiplayer and commands.
"""
from narupa.core import NarupaServer
from narupa.essd import DiscoveryServer, ServiceHub
from narupa.multiplayer.multiplayer_service import MultiplayerService, MULTIPLAYER_SERVICE_NAME
from narupa.protocol.multiplayer import add_MultiplayerServicer_to_server


class NarupaApplicationServer:
    """
    Provides a Narupa server for typical applications, with local area network discovery provided by
    ESSD, multiplayer configuration and a command service.

    """

    def __init__(self, name="Narupa Server", address=None, port=0):
        self._server = NarupaServer(address=address, port=port)
        self._discovery = DiscoveryServer()
        self._service_hub = ServiceHub(name=name,
                                       address=self._server.address,
                                       port=self._server.port)
        self._services = set()
        self._setup_multiplayer()

    @property
    def name(self) -> str:
        """
        Name of the server.
        :return: The name of the server.
        """
        return self._service_hub.name

    @property
    def address(self) -> str:
        """
        Address of the server.
        :return: Address of the server.
        """
        return self._server.address

    @property
    def port(self) -> int:
        """
        Server port.
        :return: Port of the server.
        """
        return self._server.port

    def close(self):
        """
        Close the application server and all services.
        """
        self._discovery.close()
        for service in self._services:
            service.close()
        self._server.close()

    def add_service(self, service_name, service, service_registration_method):
        """
        Adds a gRPC service to the server and discovery.
        :param service_name: Name of the service.
        :param service: Service implementation
        :param service_registration_method: Method to register service.
        """
        # TODO this seems a bit low level, but we need a way to add services to both the server and discovery.
        # TODO package up a service, service name and method?
        self._server.add_service(service, add_service_method=service_registration_method)
        self._service_hub.add_service(service_name, self._server.port)
        self._update_discovery_services()

    def _setup_multiplayer(self):
        self._multiplayer = MultiplayerService()
        self.add_service(MULTIPLAYER_SERVICE_NAME, self._multiplayer, add_MultiplayerServicer_to_server)

    def _update_discovery_services(self):
        self._discovery.unregister_service(self._service_hub)
        self._discovery.register_service(self._service_hub)
