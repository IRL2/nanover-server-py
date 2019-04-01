import time
from concurrent import futures
from queue import Queue

import grpc

from narupa.protocol.instance.molecule_provider_pb2_grpc import MoleculeProviderServicer, MoleculeProviderStub, add_MoleculeProviderServicer_to_server

from narupa.protocol.delimiter_pb2 import START, END

from narupa.protocol.instance.get_topology_pb2 import GetTopologyResponse, GetTopologyRequest

from narupa.protocol.instance.get_frame_pb2 import GetFrameResponse, GetFrameRequest

from narupa.protocol.trajectory.frame_pb2 import FrameData

from narupa.protocol.topology.topology_pb2 import TopologyData

from typing import List, Callable


class NarupaTopology(object):
    packets: List[GetTopologyResponse]

    def __init__(self, frame_index: int, topology: TopologyData):
        self.packets = [
            GetTopologyResponse(frame_index=frame_index, delimiter=START),
            GetTopologyResponse(frame_index=frame_index, topology=topology),
            GetTopologyResponse(frame_index=frame_index, delimiter=END)
        ]

class NarupaFrame(object):
    packets: List[GetFrameResponse]

    def __init__(self, frame_index: int, frame: FrameData):
        self.packets = [
            GetFrameResponse(frame_index=frame_index, delimiter=START),
            GetFrameResponse(frame_index=frame_index, frame=frame),
            GetFrameResponse(frame_index=frame_index, delimiter=END)
        ]


class NarupaClient(MoleculeProviderStub):

    def __init__(self):
        self.channel = grpc.insecure_channel('localhost:50051')
        MoleculeProviderStub.__init__(self, self.channel)
        self.executor = futures.ThreadPoolExecutor(max_workers=10)

    def topology_stream(self, callback: Callable[[List[GetTopologyResponse]], None]):
        response_list: List[GetTopologyResponse] = None
        for reply in self.SubscribeTopology(GetTopologyRequest()):
            if reply.delimiter == START:
                response_list = []
            elif reply.delimiter == END:
                callback(response_list)
                response_list = None
            else:
                response_list.append(reply)

    def frame_stream(self, callback: Callable[[List[GetFrameResponse]], None]):
        response_list: List[GetFrameResponse]
        for reply in self.SubscribeFrame(GetFrameResponse()):
            if reply.delimiter == START:
                response_list = []
            elif reply.delimiter == END:
                callback(response_list)
                response_list = None
            else:
                response_list.append(reply)

    def subscribe_topology(self, callback: Callable[[List[GetTopologyResponse]], None]):
        self.executor.submit(self.topology_stream, callback)

    def subscribe_frames(self, callback: Callable[[List[GetFrameResponse]], None]):
        self.executor.submit(self.frame_stream, callback)


class NarupaServer(object):

    def __init__(self):
        self.server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))
        self.instance_service = NarupaInstanceService()
        add_MoleculeProviderServicer_to_server(
            self.instance_service, self.server)
        self.server.add_insecure_port('[::]:50051')
        self.server.start()

    def send_topology(self, frame_index: int, topology_data: TopologyData):
        self.instance_service.send_topology(NarupaTopology(frame_index, topology_data))

    def send_frame(self, frame_index: int, frame_data : FrameData):
        self.instance_service.send_frame(NarupaFrame(frame_index, frame_data))


if __name__ == "__main__":

    server = NarupaServer()

    client = NarupaClient()

    def on_topology(list):
        print("Topology")
        print(list)

    def on_frame(list):
        print("Frame")
        print(list)

    client.subscribe_topology(on_topology)
    client.subscribe_frames(on_frame)

    time.sleep(1)

    i = 0

    while True:
        print("Loop")
        server.send_topology(i, None)
        server.send_frame(i, None)
        i += 1
        time.sleep(1.0 / 30.0)
