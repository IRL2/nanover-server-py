from concurrent import futures
from typing import Optional
from concurrent.futures import Future

import grpc

from narupa.core import get_requested_port_or_default, DEFAULT_CONNECT_ADDRESS
from narupa.protocol.trajectory import TrajectoryServiceStub, GetFrameRequest
from narupa.trajectory import FrameData
from narupa.trajectory.frame_server import DEFAULT_PORT


class FrameClient:

    def __init__(self, *, address:Optional[str]=None, port: Optional[int]=None):
        if address is None:
            address = DEFAULT_CONNECT_ADDRESS
        port = get_requested_port_or_default(port, DEFAULT_PORT)
        self.channel = grpc.insecure_channel("{0}:{1}".format(address, port))
        self.stub = TrajectoryServiceStub(self.channel)
        self.threads = futures.ThreadPoolExecutor(max_workers=10)

    def subscribe_frames_async(self, callback) -> Future:
        return self.threads.submit(self.subscribe_frames_blocking, callback)

    def subscribe_frames_blocking(self, callback):
        for response in self.stub.SubscribeFrames(GetFrameRequest()):
            try:
                callback(frame_index=response.frame_index, frame=FrameData(response.frame))
            except Exception as e:
                print(e)

    def subscribe_last_frames_async(self, callback) -> Future:
        return self.threads.submit(self.subscribe_frames_blocking, callback)

    def subscribe_last_frames_blocking(self, callback):
        for response in self.stub.SubscribeLatestFrames(GetFrameRequest()):
            callback(frame_index=response.frame_index, frame=FrameData(response.frame))

    def close(self):
        self.channel.close()
        self.threads.shutdown(wait=False)
