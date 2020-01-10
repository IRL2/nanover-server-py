from narupa.app import NarupaApplicationServer
from narupa.core import NarupaServer
from narupa.essd import DiscoveryServer
from narupa.protocol.trajectory import add_TrajectoryServiceServicer_to_server
from narupa.trajectory import FramePublisher, FRAME_SERVICE_NAME, FrameData


class NarupaFrameServer(NarupaApplicationServer):
    """

    Application-level class for implementing a Narupa frame server, something that publishes frames
    :class:`FrameData` that can be consumed, e.g. simulation trajectories.

    Example
    =======

    >>> with NarupaFrameServer.basic_server() as frame_server:
    ...     frame_publisher = frame_server.frame_publisher
    ...     example_frame = FrameData() # A simple frame representing two particles.
    ...     example_frame.particle_positions = [0,0,0,1,1,1]
    ...     example_frame.particle_count = 2
    ...     frame_publisher.send_frame(frame_index=0, frame_data=example_frame)

    """
    def __init__(self, server: NarupaServer, discovery: DiscoveryServer, name="Narupa Frame Server"):
        super().__init__(server, discovery, name)
        self._setup_frame_publisher()

    @property
    def frame_publisher(self) -> FramePublisher:
        """
        The frame publisher attached to this server. Use it to publish frames for consumption by
        Narupa frame clients.

        :return: The :class:`FramePublisher` attached to this server.
        """
        return self._frame_publisher

    def _setup_frame_publisher(self):
        self._frame_publisher = FramePublisher()
        self.add_service(FRAME_SERVICE_NAME, self._frame_publisher, add_TrajectoryServiceServicer_to_server)
