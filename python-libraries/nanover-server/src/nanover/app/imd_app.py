"""
Module providing an implementation of an NanoVer iMD application, for publishing
simulations and trajectories for consumption by clients that can be interacted
with in real-time through biasing potentials.

"""

from ssl import SSLContext
from typing import Any

from nanover.app.app_server import NanoverApplicationServer
from nanover.essd import DiscoveryServer
from nanover.imd import ImdStateWrapper
from nanover.app.multiuser import add_multiuser_commands
from nanover.websocket.server import DEFAULT_NANOVER_PORT


class NanoverImdApplication(NanoverApplicationServer):
    """
    Application-level class for implementing a NanoVer iMD server, something that publishes
    :class:`FrameData` that can be consumed, e.g. simulation trajectories, and can receive
    interactive forces in real-time, allowing the simulation to be biased.
    """

    DEFAULT_SERVER_NAME: str = "NanoVer iMD Server"
    _imd_state: ImdStateWrapper
    _server_ws: Any | None = None

    @classmethod
    def basic_server(
        cls,
        *,
        name: str | None = None,
        address: str | None = None,
        port: int | None = None,
        ssl: SSLContext | None = None,
    ):
        """
        Initialises a basic NanoVer application server with default settings,
        with a default unencrypted server and ESSD discovery server for
        finding it on a local area network.

        :param name: Name of the server for the purposes of discovery.
        :param address: The address at which to bind the server to. If none given,
            the default address of
        :param port: Optional port on which to run the NanoVer server. If none given,
            will attempt to use the default port (38801) else a random port will be used.
        :param ssl: Optional SSLContext that if provided is used to host an additional
            secure websocket server on the WSS protocol.
        :return: An instantiation of a basic NanoVer server, registered with an
            ESSD discovery server.
        """
        port = DEFAULT_NANOVER_PORT if port is None else port

        app_server = super().basic_server(name=name, address=address)
        app_server.serve_websocket(ssl=ssl, port=port)
        return app_server

    def __init__(
        self,
        *,
        discovery: DiscoveryServer | None = None,
        name: str | None = None,
        address: str | None = None,
    ):
        super().__init__(discovery=discovery, name=name, address=address)
        self._imd_state = ImdStateWrapper(self._state_service.state_dictionary)
        add_multiuser_commands(self)

    @property
    def port(self) -> int | None:
        """
        Returns first available Websocket port. 
        Insecure port will be returned if it exists else reports same as `secure_port`.
        """
        return self._server_ws.ws_port if self._server_ws is not None else self._server_ws.wss_port

    @property
    def secure_port(self) -> int | None:
        """Returns SSL wrapped Websocket port, if available."""
        return self._server_ws.wss_port if self._server_ws is not None else None

    def close(self):
        if self._server_ws is not None:
            self._server_ws.close()
        super().close()

    def serve_websocket(self, *, insecure=True, ssl: SSLContext | None = None, port: int = 0):
        from nanover.websocket.server import WebSocketServer

        self._server_ws = WebSocketServer.basic_server(self, insecure=insecure, ssl=ssl, port=port)
        return self._server_ws

    @property
    def imd(self) -> ImdStateWrapper:
        """
        The iMD service attached to this application. Use it to access interactive forces sent
        by clients, so they can be applied to a simulation.

        :return: An :class:`ImdStateWrapper` for tracking interactions.
        """
        return self._imd_state
