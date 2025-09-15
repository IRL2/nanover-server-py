import traceback
from concurrent.futures import ThreadPoolExecutor
from time import perf_counter_ns

import msgpack
from websockets.sync.client import connect, ClientConnection

from nanover.recording.writing import write_entry, write_header
from nanover.state.state_service import dictionary_change_to_state_update
from nanover.trajectory.frame_data import FRAME_INDEX
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.convert import (
    convert_dict_frame_to_grpc_frame,
    unpack_dict_frame,
)
from nanover.protocol.trajectory import GetFrameResponse


def record_from_server(address, trajectory_file, state_file):
    """
    Connect to the given host:port and record trajectory frames and state updates to files

    :param address: String host:port of server to connect to
    :param trajectory_file: File to write trajectory frames to
    :param state_file: File to write state updates to
    :return:
    """

    frame_out = open(trajectory_file, "wb")
    state_out = open(state_file, "wb")

    write_header(frame_out)
    write_header(state_out)

    start_time = perf_counter_ns()

    def timestamp():
        return int((perf_counter_ns() - start_time) / 1000)

    def recv_frame(message: dict):
        message = unpack_dict_frame(message)
        response = GetFrameResponse(
            frame_index=message[FRAME_INDEX],
            frame=convert_dict_frame_to_grpc_frame(message).raw,
        )
        write_entry(frame_out, timestamp(), response)

    def recv_state(message: dict):
        change = DictionaryChange(
            updates=message["updates"],
            removals=message["removals"],
        )
        update = dictionary_change_to_state_update(change)
        write_entry(state_out, timestamp(), update)

    def listen(websocket: ClientConnection):
        with frame_out, state_out:
            try:
                for data in websocket:
                    message = msgpack.unpackb(data)
                    frame = message.get("frame", None)
                    state = message.get("state", None)

                    if frame is not None:
                        recv_frame(frame)
                    if state is not None:
                        recv_state(state)
            except Exception:
                traceback.print_exc()
                raise
            finally:
                websocket.close()

    websocket = connect(address)
    executor = ThreadPoolExecutor(max_workers=1)
    executor.submit(listen, websocket)

    return executor, websocket
