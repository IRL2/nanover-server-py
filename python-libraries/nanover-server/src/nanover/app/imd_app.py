"""
Module providing an implementation of an NanoVer iMD application, for publishing
simulations and trajectories for consumption by clients that can be interacted
with in real-time through biasing potentials.

"""

from ssl import SSLContext
from typing import Optional, Any

from nanover.app.frame_app import NanoverFrameApplication
from nanover.core import NanoverServer
from nanover.essd import DiscoveryServer
from nanover.imd import ImdStateWrapper, IMD_SERVICE_NAME


class NanoverImdApplication(NanoverFrameApplication):
    """
    Application-level class for implementing a NanoVer iMD server, something that publishes
    :class:`FrameData` that can be consumed, e.g. simulation trajectories, and can received
    interactive forces in real-time, allowing the simulation to be biased.

    >>> with NanoverImdApplication.basic_server() as app: # fire up interactive molecular dynamics
    ...     with NanoverImdClient() as client:
    ...         client.interactions # print any active interactions (in this case, none).
    {}

    """

    DEFAULT_SERVER_NAME: str = "NanoVer iMD Server"
    _imd_state: ImdStateWrapper
    _server_ws: Optional[Any] = None

    @classmethod
    def basic_server(
        cls,
        name: Optional[str] = None,
        address: Optional[str] = None,
        port: Optional[int] = None,
        *,
        ssl: Optional[SSLContext] = None,
    ):
        """
        Initialises a basic NanoVer application server with default settings,
        with a default unencrypted server and ESSD discovery server for
        finding it on a local area network.

        :param name: Name of the server for the purposes of discovery.
        :param address: The address at which to bind the server to. If none given,
            the default address of
        :param port: Optional port on which to run the NanoVer server. If none given,
            default port will be used.
        :return: An instantiation of a basic NanoVer server, registered with an
            ESSD discovery server.
        """
        app_server = super().basic_server(name=name, address=address, port=port)
        app_server.serve_websocket(ssl=ssl)
        return app_server

    def __init__(
        self,
        server: NanoverServer,
        discovery: Optional[DiscoveryServer] = None,
        name: Optional[str] = None,
    ):
        super().__init__(server, discovery, name)
        self._setup_imd()

    def close(self):
        if self._server_ws is not None:
            self._server_ws.close()
        super().close()

    def serve_websocket(self, *, insecure=True, ssl: Optional[SSLContext] = None):
        from nanover.websocket.server import WebSocketServer

        self._server_ws = WebSocketServer.basic_server(self, insecure=insecure, ssl=ssl)
        return self._server_ws

    @property
    def imd(self) -> ImdStateWrapper:
        """
        The iMD service attached to this application. Use it to access interactive forces sent
        by clients, so they can be applied to a simulation.

        :return: An :class:`ImdStateWrapper` for tracking interactions.
        """
        return self._imd_state

    def _setup_imd(self):
        self._imd_state = ImdStateWrapper(self.server._state_service.state_dictionary)
        self._add_service_entry(IMD_SERVICE_NAME, self.server.port)
