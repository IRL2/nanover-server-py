"""
Module providing an out-of-the-box Narupa application server,
with an underyling gRPC server, discovery, multiplayer and commands.
"""
import getpass
from typing import Tuple, Optional

from narupa.core import NarupaServer, DEFAULT_SERVE_ADDRESS
from narupa.essd import DiscoveryServer, ServiceHub
from narupa.multiplayer.multiplayer_service import MultiplayerService, MULTIPLAYER_SERVICE_NAME
from narupa.protocol.multiplayer import add_MultiplayerServicer_to_server

DEFAULT_NARUPA_PORT = 38801
DEFAULT_SERVER_NAME = "Narupa Server"


def start_default_server_and_discovery(
        address: Optional[str] = None,
        port: Optional[int] = None) -> Tuple[NarupaServer, DiscoveryServer]:
    """
    Utility method for creating a default Narupa server along with ESSD discovery.

    :param: address: Address to run the server at. If nothing is passed, the default
        address of all interfaces will be used.
    :param: port: Port to run the server on, if nothing is passed, the default
        Narupa port will be used. The value of zero should be passed to let the OS
        pick a free port.
    :return: tuple of Narupa server and ESSD discovery.
    """
    address = address or DEFAULT_SERVE_ADDRESS
    if port is None:
        port = DEFAULT_NARUPA_PORT
    try:
        server = NarupaServer(address=address, port=port)
    except IOError:
        if port == DEFAULT_NARUPA_PORT:
            raise IOError(f'Could not start a server at the default port ({port}). Is another Narupa server running? '
                          f'Use port=0 to let the OS find a free port')
        raise
    discovery = DiscoveryServer()
    return server, discovery


class NarupaApplicationServer:
    """
    Provides a convenient Narupa server for typical applications, with local
    area network discovery provided by ESSD, multiplayer configuration and a
    command service.

    Use this a base for building specific applications by inheriting from it
    and attaching additional services.
    """

    def __init__(
            self,
            server: NarupaServer,
            discovery: Optional[DiscoveryServer] = None,
            name: Optional[str] = None,
    ):
        if name is None:
            name = qualified_server_name(DEFAULT_SERVER_NAME)
        self._server = server
        self._discovery = discovery
        self._service_hub = ServiceHub(name=name,
                                       address=self._server.address,
                                       port=self._server.port)
        self._services = set()
        self._setup_multiplayer()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @classmethod
    def basic_server(
            cls,
            name: Optional[str] = None,
            address: Optional[str] = None,
            port: Optional[int] = None,
    ):
        """
        Initialises a basic Narupa application server with default settings,
        with a default unencrypted server and ESSD discovery server for
        finding it on a local area network.

        :param name: Name of the server for the purposes of discovery.
        :param address: The address at which to bind the server to. If none given,
            the default address of
        :param port: Optional port on which to run the Narupa server. If none given,
            default port will be used.
        :return: An instantiation of a basic Narupa server, registered with an
            ESSD discovery server.
        """
        server, discovery = start_default_server_and_discovery(address=address, port=port)
        return cls(server, discovery, name)

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

    @property
    def server(self) -> NarupaServer:
        """
        The underlying Narupa server for this application.
        One can use this to manage commands and services.
        :return: The Narupa server.
        """
        # TODO expose command api directly?
        return self._server

    @property
    def running_discovery(self) -> bool:
        """
        Indicates whether a discovery service is running or not.
        :return: True if discovery is available, False otherwise.
        """
        return self.discovery is not None

    @property
    def discovery(self) -> Optional[DiscoveryServer]:
        """
        The discovery service that can be used to allow clients to find services hosted by this application.
        :return: The discovery service, or None if no discovery has been set up.

        Services added directly to the server running on this application via :fun:`NarupaApplicationServer.add_service`
        are automatically added to this discovery service.

        Accessing the discovery service directly enables one to register their own server that may be running
        separately to the core application.
        """
        return self._discovery

    def close(self):
        """
        Close the application server and all services.
        """
        if self.running_discovery:
            self._discovery.close()
        for service in self._services:
            service.close()
        self._server.close()

    def add_service(self, service):
        """
        Adds a gRPC service to the server and broadcast it on discovery.
        :param service: Service implementation
        """
        self._server.add_service(service)
        self._service_hub.add_service(service.name, self._server.port)
        self._services.add(service)
        if self.running_discovery:
            self._update_discovery_services()

    def _setup_multiplayer(self):
        self._multiplayer = MultiplayerService()
        self.add_service(self._multiplayer)

    def _update_discovery_services(self):
        try:
            self._discovery.unregister_service(self._service_hub)
        except KeyError:
            pass
        self._discovery.register_service(self._service_hub)


def qualified_server_name(base_name: str):
    """
    Prefixes the given server name with identifying information of the machine
    running it.
    """
    username = getpass.getuser()  # OS agnostic method that uses a few different metrics to get the username
    return f'{username}: {base_name}'

