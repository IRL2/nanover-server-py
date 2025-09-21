from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, asdict, field
from time import perf_counter_ns
from typing import IO, Iterable, Any
from zipfile import ZipFile

import msgpack
from websockets import ConnectionClosed
from websockets.sync.client import connect, ClientConnection

from nanover.recording.utilities import iter_recording_max
from nanover.state.state_service import state_update_to_dictionary_change
from nanover.trajectory.frame_data import FRAME_INDEX, FrameData
from nanover.websocket.client.base_client import MAX_MESSAGE_SIZE
from nanover.websocket.convert import convert_grpc_frame_to_dict_frame, pack_dict_frame, pack_grpc_frame

RECORDING_INDEX_FILENAME = "index.msgpack"
RECORDING_MESSAGES_FILENAME = "messages.msgpack"


@dataclass(kw_only=True)
class MessageEvent:
    timestamp: int
    data: bytes


@dataclass(kw_only=True)
class RecordingIndexEntry:
    offset: int
    length: int
    timestamp: int
    metadata: dict[str, Any] = field(default_factory=dict)

    def read_from(self, io: IO[bytes]):
        io.seek(self.offset)
        return io.read(self.length)


def read_recording_v2(in_path):
    with ZipFile(in_path) as archive:
        with archive.open(RECORDING_INDEX_FILENAME) as index_file:
            index = msgpack.unpackb(index_file.read())
        with archive.open(RECORDING_MESSAGES_FILENAME) as received_file:
            print(index[50])

            entry = RecordingIndexEntry(**index[50])
            data = entry.read_from(received_file)

            message = msgpack.unpackb(data)
            print(message)


def convert_old_recording(out_path, *, traj, state):
    record_messages(out_path, message_events_from_old_recording(traj=traj, state=state))


def record_from_server(address, out_path):
    websocket = connect(address, max_size=MAX_MESSAGE_SIZE)
    executor = ThreadPoolExecutor(max_workers=1)
    executor.submit(record_messages, out_path, message_events_from_websocket(websocket))

    return executor, websocket


def message_events_from_websocket(websocket: ClientConnection):
    start_time = perf_counter_ns()

    def get_timestamp():
        return int((perf_counter_ns() - start_time) / 1000)

    try:
        for data in websocket:
            if isinstance(data, bytes):
                yield MessageEvent(timestamp=get_timestamp(), data=data)
    except ConnectionClosed:
        pass
    finally:
        websocket.close()


def message_events_from_old_recording(*, traj, state):
    for event in iter_recording_max(traj=traj, state=state):
        if event.next_frame_event is not None:
            frame_response = event.next_frame_event.message
            frame = pack_grpc_frame(FrameData(frame_response.frame))
            frame[FRAME_INDEX] = frame_response.frame_index
            message = {"frame": frame}
            data = msgpack.packb(message)
            yield MessageEvent(timestamp=event.timestamp, data=data)
        if event.next_state_event is not None:
            change = state_update_to_dictionary_change(event.next_state_event.message)
            message = {
                "state": {"updates": change.updates, "removals": list(change.removals)}
            }
            data = msgpack.packb(message)
            yield MessageEvent(timestamp=event.timestamp, data=data)


def record_messages(out_path: str, messages: Iterable[MessageEvent]):
    with ZipFile(out_path, "w") as archive:
        entries: list[RecordingIndexEntry] = []
        try:
            with archive.open(RECORDING_MESSAGES_FILENAME, "w") as messages_file:
                offset = 0
                for event in messages:
                    entry = RecordingIndexEntry(
                        offset=offset,
                        length=len(event.data),
                        timestamp=event.timestamp,
                    )
                    entries.append(entry)
                    messages_file.write(event.data)
                    offset += entry.length
        finally:
            with archive.open(RECORDING_INDEX_FILENAME, "w") as index_file:
                data = [asdict(entry) for entry in entries]
                index_file.write(msgpack.packb(data))
