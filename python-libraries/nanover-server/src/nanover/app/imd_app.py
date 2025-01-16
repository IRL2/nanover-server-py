"""
Module providing an implementation of an NanoVer iMD application, for publishing
simulations and trajectories for consumption by clients that can be interacted
with in real-time through biasing potentials.

"""

from typing import Optional

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

    def __init__(
        self,
        server: NanoverServer,
        discovery: Optional[DiscoveryServer] = None,
        name: Optional[str] = None,
    ):
        super().__init__(server, discovery, name)
        self._setup_imd()

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
