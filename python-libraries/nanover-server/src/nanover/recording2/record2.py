from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, asdict, field
from os import PathLike
from time import perf_counter_ns
from typing import IO, Iterable, Any, BinaryIO
from zipfile import ZipFile

import msgpack
from websockets import ConnectionClosed
from websockets.sync.client import connect, ClientConnection

from nanover.trajectory.frame_data import FRAME_INDEX
from nanover.websocket.client.base_client import MAX_MESSAGE_SIZE

RECORDING_INDEX_FILENAME = "index.msgpack"
RECORDING_MESSAGES_FILENAME = "messages.msgpack"


@dataclass(kw_only=True)
class MessageEvent:
    timestamp: int
    message: dict[str, Any]


@dataclass(kw_only=True)
class RecordingIndexEntry:
    offset: int
    length: int
    metadata: dict[str, Any] = field(default_factory=dict)

    def read_from(self, io: IO[bytes]):
        io.seek(self.offset)
        return io.read(self.length)


class MessageZipReader:
    @classmethod
    def from_path(cls, path: PathLike[str]):
        """
        Read a recording from a filepath.
        """
        return cls(ZipFile(path, "r"))

    @classmethod
    def from_io(cls, io: BinaryIO):
        """
        Read a recording from binary data source.
        """
        return cls(ZipFile(io, "r"))

    def __init__(self, zipfile: ZipFile):
        self.zipfile = zipfile
        self.messagesfile = zipfile.open(RECORDING_MESSAGES_FILENAME)
        self.index = parse_index(zipfile)

    def close(self):
        self.zipfile.close()

    def get_message_at_index(self, index):
        return self.get_message_from_entry(self.index[index])

    def get_message_from_entry(self, entry: RecordingIndexEntry):
        return msgpack.unpackb(entry.read_from(self.messagesfile))

    def iter_messages(self):
        for entry in self:
            yield entry.timestamp, self.get_message_from_entry(entry)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        yield from self.index

    def __len__(self):
        return len(self.index)

    def __getitem__(self, index):
        return self.index[index]


class NanoVerRecordingReader(MessageZipReader):
    def iter_frames(self):
        for entry in self.index:
            if "frame" in entry.metadata["types"]:
                message = self.get_message_from_entry(entry)
                yield message["frame"]


def parse_index(zipfile: ZipFile):
    with zipfile.open(RECORDING_INDEX_FILENAME) as index_file:
        return [
            RecordingIndexEntry(**entry) for entry in msgpack.unpackb(index_file.read())
        ]


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
                yield MessageEvent(
                    timestamp=get_timestamp(),
                    message=msgpack.unpackb(data),
                )
    except ConnectionClosed:
        pass
    finally:
        websocket.close()


def record_messages(out_path: str, messages: Iterable[MessageEvent]):
    with ZipFile(out_path, "w") as archive:
        entries: list[RecordingIndexEntry] = []
        try:
            with archive.open(
                RECORDING_MESSAGES_FILENAME, "w", force_zip64=True
            ) as messages_file:
                offset = 0
                for event in messages:
                    data = msgpack.packb(event.message)
                    metadata = generate_metadata(event.message)
                    metadata["timestamp"] = event.timestamp
                    entry = RecordingIndexEntry(
                        offset=offset,
                        length=len(data),
                        metadata=metadata,
                    )
                    entries.append(entry)
                    messages_file.write(data)
                    offset += entry.length
        finally:
            with archive.open(RECORDING_INDEX_FILENAME, "w") as index_file:
                data = [asdict(entry) for entry in entries]
                index_file.write(msgpack.packb(data))


def generate_metadata(message: dict[str, Any]) -> dict[str, Any]:
    metadata = {
        "types": list(message.keys()),
    }

    if "frame" in message and FRAME_INDEX in message["frame"]:
        metadata[FRAME_INDEX] = message["frame"][FRAME_INDEX]

    return metadata
