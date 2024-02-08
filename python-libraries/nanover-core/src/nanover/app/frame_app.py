"""
Module providing an implementation of an NanoVer frame-serving application, for publishing
simulations and trajectories for consumption by clients.

"""

from typing import Optional

from nanover.app import NanoverApplicationServer
from nanover.core import NanoverServer
from nanover.essd import DiscoveryServer
from nanover.trajectory import FramePublisher


class NanoverFrameApplication(NanoverApplicationServer):
    """

    Application-level class for implementing a NanoVer frame server, something that publishes
    :class:`FrameData` that can be consumed, e.g. simulation trajectories.

    Example
    =======

    >>> with NanoverFrameApplication.basic_server() as app:
    ...     frame_publisher = app.frame_publisher
    ...     example_frame = FrameData() # A simple frame representing two particles.
    ...     example_frame.particle_positions = [[0,0,0],[1,1,1]]
    ...     example_frame.particle_count = 2
    ...     frame_publisher.send_frame(0, example_frame)

    """

    DEFAULT_SERVER_NAME: str = "NanoVer Frame Server"

    def __init__(
        self,
        server: NanoverServer,
        discovery: Optional[DiscoveryServer] = None,
        name: Optional[str] = None,
    ):
        super().__init__(server, discovery, name)
        self._setup_frame_publisher()

    def close(self):
        self._frame_publisher.close()
        super().close()

    @property
    def frame_publisher(self) -> FramePublisher:
        """
        The frame publisher attached to this application. Use it to publish
        frames for consumption by NanoVer frame clients.

        :return: The :class:`FramePublisher` attached to this application.
        """
        # TODO could just expose send frame here.
        return self._frame_publisher

    def _setup_frame_publisher(self):
        self._frame_publisher = FramePublisher()
        self.add_service(self._frame_publisher)
