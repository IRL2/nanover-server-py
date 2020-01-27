"""
Module providing an implementation of an Narupa frame-serving application, for publishing
simulations and trajectories for consumption by clients.

"""
from typing import Optional

from narupa.app import NarupaApplicationServer
from narupa.app.app_server import qualified_server_name
from narupa.core import NarupaServer
from narupa.essd import DiscoveryServer
from narupa.trajectory import FramePublisher, FrameData


class NarupaFrameApplication(NarupaApplicationServer):
    """

    Application-level class for implementing a Narupa frame server, something that publishes
    :class:`FrameData` that can be consumed, e.g. simulation trajectories.

    Example
    =======

    >>> with NarupaFrameApplication.basic_server() as app:
    ...     frame_publisher = app.frame_publisher
    ...     example_frame = FrameData() # A simple frame representing two particles.
    ...     example_frame.particle_positions = [[0,0,0],[1,1,1]]
    ...     example_frame.particle_count = 2
    ...     frame_publisher.send_frame(0, example_frame)

    """
    DEFAULT_SERVER_NAME: str = "Narupa Frame Server"

    def __init__(self, server: NarupaServer,
                 discovery: Optional[DiscoveryServer] = None,
                 name: Optional[str] = None):
        super().__init__(server, discovery, name)
        self._setup_frame_publisher()

    def close(self):
        self._frame_publisher.close()
        super().close()

    @property
    def frame_publisher(self) -> FramePublisher:
        """
        The frame publisher attached to this application. Use it to publish
        frames for consumption by Narupa frame clients.

        :return: The :class:`FramePublisher` attached to this application.
        """
        # TODO could just expose send frame here.
        return self._frame_publisher

    def _setup_frame_publisher(self):
        self._frame_publisher = FramePublisher()
        self.add_service(self._frame_publisher)
