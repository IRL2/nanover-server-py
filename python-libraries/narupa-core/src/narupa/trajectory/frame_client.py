# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from concurrent.futures import Future

import grpc
from narupa.core import NarupaStubClient
from narupa.protocol.trajectory import TrajectoryServiceStub, GetFrameRequest
from narupa.trajectory import FrameData


class FrameClient(NarupaStubClient):
    
    def __init__(self, *,
                 channel: grpc.Channel,
                 make_channel_owner: bool = False):
        super().__init__(channel=channel, stub=TrajectoryServiceStub, make_channel_owner=make_channel_owner)

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
