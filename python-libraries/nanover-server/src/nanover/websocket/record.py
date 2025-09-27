import traceback
from concurrent.futures import ThreadPoolExecutor
from time import perf_counter_ns

import msgpack
from websockets.sync.client import connect, ClientConnection

from nanover.omni import OmniRunner
from nanover.recording.writing import write_entry, write_header
from nanover.state.state_service import dictionary_change_to_state_update
from nanover.trajectory import FrameData2
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.client.app_client import get_websocket_address_from_app_server
from nanover.websocket.client.base_client import MAX_MESSAGE_SIZE
from nanover.trajectory.convert import (
    unpack_dict_frame,
    convert_framedata2_to_GetFrameResponse,
)


def record_from_runner(runner: OmniRunner, trajectory_file, state_file):
    """
    Connect to the given runner and record trajectory frames and state updates to files

    :param runner: OmniRunner instance to connect to
    :param trajectory_file: File to write trajectory frames to
    :param state_file: File to write state updates to
    :return:
    """
    return record_from_server(
        address=get_websocket_address_from_app_server(runner.app_server),
        trajectory_file=trajectory_file,
        state_file=state_file,
    )


def record_from_server(address, trajectory_file, state_file):
    """
    Connect to the given address and record trajectory frames and state updates to files

    :param address: String protocol://host:port of server to connect to
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
        frame = FrameData2(unpack_dict_frame(message))
        response = convert_framedata2_to_GetFrameResponse(frame)
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

    websocket = connect(address, max_size=MAX_MESSAGE_SIZE)
    executor = ThreadPoolExecutor(max_workers=1)
    executor.submit(listen, websocket)

    return executor, websocket
