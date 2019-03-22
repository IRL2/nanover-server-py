from narupa.protocol.instance.get_frame_pb2 import GetFrameResponse
from narupa.protocol.instance.get_topology_pb2 import GetTopologyResponse
from typing import List
from narupa.protocol.trajectory.frame_pb2 import FrameData
from narupa.protocol.topology.topology_pb2 import TopologyData
from narupa.protocol.delimiter_pb2 import START, END
from queue import Queue
from narupa.protocol.instance.molecule_provider_pb2_grpc import MoleculeProviderServicer

class GetTopologyResponseCollection(object):
    """
    Represents a collection of one or more GetTopologyResponse messages. They are collected
    such that they are initialised once and sent to every client
    """

    packets: List[GetTopologyResponse]

    def __init__(self, frame_index: int, topology: TopologyData):
        self.packets = [
            GetTopologyResponse(frame_index=frame_index, topology=topology)
        ]

class GetFrameResponseCollection(object):
    packets: List[GetFrameResponse]

    def __init__(self, frame_index: int, frame: FrameData):
        self.packets = [
            GetFrameResponse(frame_index=frame_index, frame=frame)
        ]

class FramePublisher(MoleculeProviderServicer):
    """
    An implementation of a molecule provider service. Call send_frame and send_topology
    to send data to clients when called by other python code.
    """

    topology_queues: List[Queue]
    frame_queues: List[Queue]

    last_topology: GetTopologyResponseCollection
    last_frame: GetFrameResponseCollection

    def __init__(self):
        self.topology_queues = []
        self.frame_queues = []
        self.last_frame = None
        self.last_topology = None

    def SubscribeTopology(self, request, context):

        if self.last_topology is not None:
            for packet in self.last_topology.packets:
                yield packet

        queue = Queue()
        self.topology_queues.append(queue)

        while True:
            item = queue.get(True)
            for packet in item.packets:
                yield packet

    def SubscribeFrame(self, request, context):

        if self.last_frame is not None:
            for packet in self.last_frame.packets:
                yield packet

        queue = Queue()
        self.frame_queues.append(queue)

        while True:
            item = queue.get(True)
            for packet in item.packets:
                yield packet

    def send_topology(self, frame_index : int, topology: TopologyData):
        message = GetTopologyResponseCollection(frame_index, topology)
        for queue in self.topology_queues:
            queue.put(message)
        self.last_topology = message

    def send_frame(self, frame_index : int, frame: FrameData):
        message = GetFrameResponseCollection(frame_index, frame)
        for queue in self.frame_queues:
            queue.put(message)
        self.last_frame = message

