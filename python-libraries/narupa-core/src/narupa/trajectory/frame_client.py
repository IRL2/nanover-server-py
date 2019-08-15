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

    def subscribe_frames_async(self, callback, frame_interval=0) -> Future:
        return self.threads.submit(self.subscribe_frames_blocking,
                                   callback,
                                   frame_interval)

    def subscribe_frames_blocking(self, callback, frame_interval=0):
        for frame_index, frame in self.subscribe_frames_iterate(frame_interval):
            callback(frame_index=frame_index, frame=frame)

    def subscribe_frames_iterate(self, frame_interval=0):
        request = GetFrameRequest(frame_interval=frame_interval)
        for response in self.stub.SubscribeFrames(request):
            yield response.frame_index, FrameData(response.frame)

    def subscribe_last_frames_async(self, callback, frame_interval=0) -> Future:
        return self.threads.submit(self.subscribe_last_frames_blocking,
                                   callback,
                                   frame_interval)

    def subscribe_last_frames_blocking(self, callback, frame_interval=0):
        request = GetFrameRequest(frame_interval=frame_interval)
        for response in self.stub.SubscribeLatestFrames(request):
            callback(frame_index=response.frame_index,
                     frame=FrameData(response.frame))

    def close(self):
        self.channel.close()
        self.threads.shutdown(wait=False)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
