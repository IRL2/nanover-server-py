"""
Module providing an implementation of an Narupa frame-serving application, for publishing
simulations and trajectories for consumption by clients.

"""
from narupa.app import NarupaApplicationServer
from narupa.core import NarupaServer
from narupa.essd import DiscoveryServer
from narupa.protocol.trajectory import add_TrajectoryServiceServicer_to_server
from narupa.trajectory import FramePublisher, FRAME_SERVICE_NAME, FrameData


class NarupaFrameApplication(NarupaApplicationServer):
    """

    Application-level class for implementing a Narupa frame server, something that publishes
    :class:`FrameData` that can be consumed, e.g. simulation trajectories.

    Example
    =======

    >>> with NarupaFrameApplication.basic_server() as app:
    ...     frame_publisher = app.frame_publisher
    ...     example_frame = FrameData() # A simple frame representing two particles.
    ...     example_frame.particle_positions = [0,0,0,1,1,1]
    ...     example_frame.particle_count = 2
    ...     frame_publisher.send_frame(frame_index=0, frame_data=example_frame)

    """
    def __init__(self, server: NarupaServer, discovery: DiscoveryServer, name="Narupa Frame Server"):
        super().__init__(server, discovery, name)
        self._setup_frame_publisher()

    def close(self):
        self._frame_publisher.close()
        super().close()

    @property
    def frame_publisher(self) -> FramePublisher:
        """
        The frame publisher attached to this application. Use it to publish frames for consumption by
        Narupa frame clients.

        :return: The :class:`FramePublisher` attached to this application.
        """
        # TODO could just expose send frame here.
        return self._frame_publisher

    def _setup_frame_publisher(self):
        self._frame_publisher = FramePublisher()
        self.add_service(FRAME_SERVICE_NAME, self._frame_publisher, add_TrajectoryServiceServicer_to_server)
