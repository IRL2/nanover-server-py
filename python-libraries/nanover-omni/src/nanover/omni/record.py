import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import grpc

from nanover.protocol.state import StateStub, SubscribeStateUpdatesRequest
from nanover.protocol.trajectory import TrajectoryServiceStub, GetFrameRequest
from nanover.recording.writing import record_messages


from threading import Lock


def record_from_server(address, trajectory_file, state_file):
    """
    Connect to the given host:port and record trajectory frames and state updates to files
    :param address: String host:port of server to connect to
    :param trajectory_file: File to write trajectory frames to
    :param state_file: File to write state updates to
    :return:
    """
    print_lock = Lock()

    def error_handler(f):
        e = f.exception()

        if e is None:
            return

        with print_lock:
            traceback.print_exc()

    channel = grpc.insecure_channel(address)

    executor = ThreadPoolExecutor(max_workers=2)
    executor.submit(record_trajectory, trajectory_file, channel).add_done_callback(
        error_handler
    )
    executor.submit(record_state, state_file, channel).add_done_callback(error_handler)
    return executor, channel


def record_trajectory(path, channel):
    path = Path(path)
    path.parent.mkdir(exist_ok=True, parents=True)

    with path.open("wb") as io:
        stub = TrajectoryServiceStub(channel)
        request = GetFrameRequest()
        stream = stub.SubscribeLatestFrames(request)
        record_messages(io, stream)


def record_state(path, channel):
    path = Path(path)
    path.parent.mkdir(exist_ok=True, parents=True)

    with path.open("wb") as io:
        stub = StateStub(channel)
        request = SubscribeStateUpdatesRequest()
        stream = stub.SubscribeStateUpdates(request)
        record_messages(io, stream)
