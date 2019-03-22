from narupy.core import GrpcServer
from narupa.protocol.trajectory.frame_pb2 import FrameData
from narupa.protocol.topology.topology_pb2 import TopologyData
from narupa.protocol.instance.molecule_provider_pb2_grpc import add_MoleculeProviderServicer_to_server
from narupy.trajectory.frame_publisher import FramePublisher


class FrameServer(GrpcServer):

    def setup_services(self):
        super().setup_services()
        self._instance_service = FramePublisher()
        add_MoleculeProviderServicer_to_server(self._instance_service, self.server)

    def send_topology(self, frame_index: int, topology_data: TopologyData):
        self._instance_service.send_topology(frame_index, topology_data)

    def send_frame(self, frame_index: int, frame_data: FrameData):
        self._instance_service.send_frame(frame_index, frame_data)
