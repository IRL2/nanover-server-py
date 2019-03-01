# Copyright (c) Mike O'Connor, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from concurrent import futures
from typing import List, Iterator
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import MultiplayerServicer
from queue import Queue

class MultiplayerService(MultiplayerServicer):
    avatar_queues: List[Queue]


    def __init__(self):
        self.avatar_queues = []


    def SubscribeToAvatars(self, request, context):
        queue = Queue()
        self.avatar_queues.append(queue)

        while True:
            item = queue.get(True)
            for packet in item.packets:
                yield packet

    def PublishAvatar(self, request_iterator : Iterator, context):

        # TODO right publish.
        pass


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

    def send_topology(self, topology: NarupaTopology):
        for queue in self.topology_queues:
            queue.put(topology)
        self.last_topology = topology

    def send_frame(self, frame: NarupaFrame):
        for queue in self.frame_queues:
            queue.put(frame)
        self.last_frame = frame

class MultiplayerServer(object):
    """
    Server providing multiplayer synchronisation.
    """

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