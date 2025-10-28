"""
Module providing an out-of-the-box NanoVer application server supporting all the standard features.
"""

import getpass
from ssl import SSLContext
from typing import Any

from nanover.essd import DiscoveryServer, ServiceHub
from nanover.trajectory import FramePublisher
from nanover.utilities.change_buffers import DictionaryChange

from nanover.core.commands import CommandService, CommandHandler
from nanover.utilities.state_dictionary import StateDictionary
from .multiuser import add_multiuser_commands
from ..core import AppServer
from ..imd import ImdStateWrapper

DEFAULT_SERVE_ADDRESS = "[::]"
DEFAULT_NANOVER_PORT = 38801


class NanoverImdApplication(AppServer):
    """
    Application-level class for implementing a NanoVer iMD server, something that publishes
    :class:`FrameData` that can be consumed, e.g. simulation trajectories, and can receive
    interactive forces in real-time, allowing the simulation to be biased.
    """

    DEFAULT_SERVER_NAME: str = "NanoVer iMD Server"
    _imd_state: ImdStateWrapper
    _server_ws: Any | None = None

    @classmethod
    def null_server(cls):
        return cls()

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
        :param port: Optional port on which to run the NanoVer server, defaulting to 38801 (the NanoVer default port).
        :param ssl: Optional SSLContext that if provided is used to host an additional
            secure websocket server on the WSS protocol.
        :return: An instantiation of a basic NanoVer server, registered with an
            ESSD discovery server.
        """
        port = DEFAULT_NANOVER_PORT if port is None else port

        app_server = cls(name=name, address=address, discovery=DiscoveryServer())
        app_server.serve_websocket(ssl=ssl, port=port)
        return app_server

    def __init__(
        self,
        *,
        discovery: DiscoveryServer | None = None,
        name: str | None = None,
        address: str | None = None,
    ):
        if name is None:
            name = qualified_server_name(self.DEFAULT_SERVER_NAME)

        if address is None:
            address = DEFAULT_SERVE_ADDRESS

        self._address = address
        self._command_service = CommandService()
        self._state_dictionary = StateDictionary()
        self._frame_publisher = FramePublisher()

        self._discovery = discovery
        self._service_hub = ServiceHub(name=name, address=address)

        self._imd_state = ImdStateWrapper(self._state_dictionary)
        add_multiuser_commands(self)

    @property
    def name(self):
        """
        Name of the server.
        :return: The name of the server.
        """
        return self._service_hub.name

    @property
    def address(self):
        """
        Address of the server.
        :return: Address of the server.
        """
        return self._address

    @property
    def discovery(self):
        """
        The discovery service that can be used to allow clients to find services hosted by this application.
        :return: The discovery service, or None if no discovery has been set up.

        Services added directly to the server running on this application via :func:`NanoverApplicationServer.add_service`
        are automatically added to this discovery service.

        Accessing the discovery service directly enables one to register their own server that may be running
        separately to the core application.
        """
        return self._discovery

    @property
    def frame_publisher(self):
        """
        The frame publisher attached to this application. Use it to publish
        frames for consumption by NanoVer frame clients.

        :return: The :class:`FramePublisher` attached to this application.
        """
        # TODO could just expose send frame here.
        return self._frame_publisher

    @property
    def imd(self) -> ImdStateWrapper:
        """
        The iMD service attached to this application. Use it to access interactive forces sent
        by clients, so they can be applied to a simulation.

        :return: An :class:`ImdStateWrapper` for tracking interactions.
        """
        return self._imd_state

    def serve_websocket(
        self, *, insecure=True, ssl: SSLContext | None = None, port: int = 0
    ):
        from nanover.websocket.server import WebSocketServer

        self._server_ws = WebSocketServer.basic_server(
            self, insecure=insecure, ssl=ssl, port=port
        )
        return self._server_ws

    def close(self):
        """
        Close the application server and all services.
        """
        self._state_dictionary.freeze()
        self._frame_publisher.close()
        if self._discovery is not None:
            self._discovery.close()
        if self._server_ws is not None:
            self._server_ws.close()

    @property
    def commands(self):
        """
        Gets the commands available on this server.

        :return: The commands, consisting of their names, callback and default parameters.
        """
        return self._command_service.commands

    def run_command(self, name: str, arguments: dict[str, Any]):
        return self._command_service.run_command(name, arguments)

    def register_command(
        self,
        name: str,
        callback: CommandHandler,
        default_arguments: dict | None = None,
    ):
        """
        Adds a command on this server.

        :param name: Name of the command to register
        :param callback: Method to be called whenever the given command name is run by a client.
        :param default_arguments: A description of the arguments of the callback and their default values.

        :raises ValueError: Raised when a command with the same name already exists.
        """
        self._command_service.register_command(name, callback, default_arguments)

    def unregister_command(self, name):
        """
        Deletes a command from this server.

        :param name: Name of the command to delete
        """
        self._command_service.unregister_command(name)

    def lock_state(self):
        """
        Context manager for reading the current state while delaying any changes
        to it.
        """
        return self._state_dictionary.lock_content()

    def copy_state(self):
        """
        Return a shallow copy of the current state.
        """
        return self._state_dictionary.copy_content()

    def update_state(self, access_token: Any, change: DictionaryChange):
        """
        Attempts an atomic update of the shared key/value store. If any key
        cannot be updates, no change will be made.
        """
        self._state_dictionary.update_state(access_token, change)

    def clear_locks(self):
        """
        Forces the release all locks on all keys in the shared key/value store.
        """
        self._state_dictionary.clear_locks()

    @property
    def state_dictionary(self):
        return self._state_dictionary

    def add_service(self, name: str, port: int):
        self._service_hub.add_service(name, port)
        if self._discovery is not None:
            self._update_discovery_services()

    @property
    def service_hub(self):
        return self._service_hub

    def _update_discovery_services(self):
        try:
            self._discovery.unregister_service(self._service_hub)  # type: ignore
        except KeyError:
            pass
        self._discovery.register_service(self._service_hub)  # type: ignore

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def qualified_server_name(base_name: str):
    """
    Prefixes the given server name with identifying information of the machine
    running it.
    """
    username = (
        getpass.getuser()
    )  # OS agnostic method that uses a few different metrics to get the username
    return f"{username}: {base_name}"
