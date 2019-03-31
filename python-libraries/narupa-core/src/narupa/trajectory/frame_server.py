<<<<<<< HEAD:project/python/narupa-protocol/narupy/trajectory/frame_server.py
from narupa.protocol.instance.molecule_provider_pb2_grpc import add_MoleculeProviderServicer_to_server
from narupa.protocol.trajectory.frame_pb2 import FrameData
from narupy.core import GrpcServer
from narupy.trajectory.frame_publisher import FramePublisher
=======
from narupa.core import GrpcServer
from narupa.protocol.trajectory.frame_pb2 import FrameData
from narupa.protocol.topology.topology_pb2 import TopologyData
from narupa.protocol.instance.molecule_provider_pb2_grpc import add_MoleculeProviderServicer_to_server
from .frame_publisher import FramePublisher
>>>>>>> b942fec690e7bfd12b38ca051fb3a8dc750436b2:python-libraries/narupa-core/src/narupa/trajectory/frame_server.py


class FrameServer(GrpcServer):
    _trajectory_service: FramePublisher

    def setup_services(self):
        super().setup_services()
        self._trajectory_service = FramePublisher()
        add_MoleculeProviderServicer_to_server(self._trajectory_service, self.server)

    def send_frame(self, frame_index: int, frame_data: FrameData):
        self._trajectory_service.send_frame(frame_index, frame_data)
