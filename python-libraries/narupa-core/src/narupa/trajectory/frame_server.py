from typing import Optional
from narupa.core import GrpcServer
from narupa.protocol.trajectory import add_TrajectoryServiceServicer_to_server
from .frame_data import FrameData
from .frame_publisher import FramePublisher

DEFAULT_ADDRESS = '[::]'
DEFAULT_PORT = 54321


class FrameServer(GrpcServer):
    _trajectory_service: FramePublisher

    def __init__(self, *, address: Optional[str]=None, port: Optional[int]=None):
        if address is None:
            address = DEFAULT_ADDRESS
        if port is None:
            port = DEFAULT_PORT
        super().__init__(address=address, port=port)
        self._frame_count = 0

    def setup_services(self):
        super().setup_services()
        self._trajectory_service = FramePublisher()
        add_TrajectoryServiceServicer_to_server(self._trajectory_service, self.server)

    def send_frame(self, frame_index: int, frame_data: FrameData):
        self._trajectory_service.send_frame(frame_index, frame_data.raw)
        self._frame_count += 1

    @property
    def frame_count(self):
        """
        Counts how many times send_frame has been called on this publisher.
        """
        return self._frame_count
