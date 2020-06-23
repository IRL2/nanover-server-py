"""
Module providing an implementation of an Narupa iMD application, for publishing
simulations and trajectories for consumption by clients that can be interacted
with in real-time through biasing potentials.

"""
from typing import Optional

from narupa.app import NarupaImdClient
from narupa.app.frame_app import NarupaFrameApplication
from narupa.core import NarupaServer
from narupa.essd import DiscoveryServer
from narupa.imd import ImdService, IMD_SERVICE_NAME


class NarupaImdApplication(NarupaFrameApplication):
    """

    Application-level class for implementing a Narupa iMD server, something that publishes
    :class:`FrameData` that can be consumed, e.g. simulation trajectories, and can received
    interactive forces in real-time, allowing the simulation to be biased.

    >>> with NarupaImdApplication.basic_server() as app: # fire up interactive molecular dynamics
    ...     with NarupaImdClient() as client:
    ...         client.interactions # print any active interactions (in this case, none).
    {}

    """
    DEFAULT_SERVER_NAME: str = "Narupa iMD Server"
    _imd_service: ImdService

    def __init__(self, server: NarupaServer,
                 discovery: Optional[DiscoveryServer] = None,
                 name: Optional[str] = None):
        super().__init__(server, discovery, name)
        self._setup_imd()

    @property
    def imd(self) -> ImdService:
        """
        The iMD service attached to this application. Use it to access interactive forces sent
        by clients, so they can be applied to a simulation.

        :return: The :class:`ImdService` attached to this application.
        """
        return self._imd_service

    def close(self):
        self._imd_service.close()
        super().close()

    def _setup_imd(self):
        self._imd_service = ImdService(self.server._state_service.state_dictionary)
        self._add_service_entry(IMD_SERVICE_NAME, self.server.port)
