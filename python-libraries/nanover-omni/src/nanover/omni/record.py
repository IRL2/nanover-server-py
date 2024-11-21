import time
from concurrent.futures import ThreadPoolExecutor

import grpc

from nanover.protocol.state import StateStub, SubscribeStateUpdatesRequest
from nanover.protocol.trajectory import TrajectoryServiceStub, GetFrameRequest
from nanover.recording.writing import record_messages


def record_from_server(address, trajectory_file, state_file):
    """
    Connect to the given host:port and record trajectory frames and state updates to files
    :param address: String host:port of server to connect to
    :param trajectory_file: File to write trajectory frames to
    :param state_file: File to write state updates to
    :return:
    """
    channel = grpc.insecure_channel(address)
    executor = ThreadPoolExecutor(max_workers=2)
    executor.submit(record_trajectory, trajectory_file, channel)
    executor.submit(record_state, state_file, channel)
    return executor, channel


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


def main():
    executor, channel = record_from_server("localhost:38801", "test.traj", "test.state")

    print("recording from server, press ctrl+c to finish")

    try:
        while True:
            time.sleep(0.01)
    except KeyboardInterrupt:
        pass

    channel.close()

    print("recording finished")


if __name__ == "__main__":
    main()
