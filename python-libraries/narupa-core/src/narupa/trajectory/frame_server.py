from narupa.core import GrpcServer
from narupa.protocol.instance import add_TrajectoryServiceServicer_to_server
from narupa.protocol.trajectory import FrameData
from narupa.trajectory import FramePublisher


class FrameServer(GrpcServer):
    _trajectory_service: FramePublisher

    def setup_services(self):
        super().setup_services()
        self._trajectory_service = FramePublisher()
        add_TrajectoryServiceServicer_to_server(self._trajectory_service, self.server)

    def send_frame(self, frame_index: int, frame_data: FrameData):
        self._trajectory_service.send_frame(frame_index, frame_data)
