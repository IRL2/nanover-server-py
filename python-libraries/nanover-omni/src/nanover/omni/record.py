from concurrent.futures import ThreadPoolExecutor

import grpc

from nanover.protocol.state import StateStub, SubscribeStateUpdatesRequest
from nanover.protocol.trajectory import TrajectoryServiceStub, GetFrameRequest
from nanover.recording.writing import record_messages


def record_from_server(address, trajectory_file, state_file):
    channel = grpc.insecure_channel(address)
    executor = ThreadPoolExecutor(max_workers=2)
    executor.submit(record_trajectory, trajectory_file, channel)
    executor.submit(record_state, state_file, channel)
    return executor


def record_trajectory(path, channel):
    with open(path, "wb") as io:
        stub = TrajectoryServiceStub(channel)
        request = GetFrameRequest()
        stream = stub.SubscribeLatestFrames(request)
        record_messages(io, stream)


def record_state(path, channel):
    with open(path, "wb") as io:
        stub = StateStub(channel)
        request = SubscribeStateUpdatesRequest()
        stream = stub.SubscribeStateUpdates(request)
        record_messages(io, stream)
