import traceback
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from time import perf_counter_ns
from zipfile import ZipFile

import msgpack
from websockets import ConnectionClosed
from websockets.sync.client import connect, ClientConnection

from nanover.omni import OmniRunner
from nanover.recording.utilities import iter_recording_max
from nanover.recording.writing import write_entry, write_header
from nanover.state.state_service import dictionary_change_to_state_update
from nanover.trajectory.frame_data import FRAME_INDEX
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.client.app_client import get_websocket_address_from_app_server
from nanover.websocket.client.base_client import MAX_MESSAGE_SIZE
from nanover.websocket.convert import (
    convert_dict_frame_to_grpc_frame,
    unpack_dict_frame, pack_dict_frame, convert_grpc_frame_to_dict_frame,
)
from nanover.protocol.trajectory import GetFrameResponse


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


@dataclass(kw_only=True)
class RecordingIndexMessageEntry:
    offset: int
    length: int
    timestamp: int


def convert_recording(out_path, *, traj, state):
    for event in iter_recording_max(traj=traj, state=state):
        if event.next_frame_event is not None:
            frame_response = event.next_frame_event.message
            frame = convert_grpc_frame_to_dict_frame(frame_response.frame)
            frame[FRAME_INDEX] = frame_response.index
            message = {"frame": pack_dict_frame(frame)}


def read_recording_v2(in_path):
    with ZipFile(in_path) as archive:
        with archive.open("index.msgpack") as index_file:
            index = msgpack.unpackb(index_file.read())
        with archive.open("received.msgpack") as received_file:
            print(index[50])

            offset = index[50]["offset"]
            length = index[51]["offset"] - offset

            received_file.seek(offset)
            data = received_file.read(length)

            message = msgpack.unpackb(data)
            print(message)


def record_from_server_v2(address, out_path):
    start_time = perf_counter_ns()

    def get_timestamp():
        return int((perf_counter_ns() - start_time) / 1000)

    def listen(websocket: ClientConnection):
        with ZipFile(out_path, "w") as archive:
            entries: list[RecordingIndexMessageEntry] = []

            with archive.open("received.msgpack", "w") as messages_file:
                offset = 0

                try:
                    for data in websocket:
                        timestamp = get_timestamp()
                        entries.append(
                            RecordingIndexMessageEntry(
                                offset=offset,
                                length=len(data),
                                timestamp=timestamp,
                            )
                        )
                        messages_file.write(data)
                        offset += len(data)
                except ConnectionClosed:
                    pass
                except Exception:
                    traceback.print_exc()
                    raise
                finally:
                    websocket.close()

            with archive.open("index.msgpack", "w") as index_file:
                data = [{
                    "offset": entry.offset,
                    "timestamp": entry.timestamp,
                    "length": entry.length,
                } for entry in entries]
                index_file.write(msgpack.packb(data))

    websocket = connect(address, max_size=MAX_MESSAGE_SIZE)
    executor = ThreadPoolExecutor(max_workers=1)
    executor.submit(listen, websocket)

    return executor, websocket


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

    websocket = connect(address, max_size=MAX_MESSAGE_SIZE)
    executor = ThreadPoolExecutor(max_workers=1)
    executor.submit(listen, websocket)

    return executor, websocket
