from dataclasses import dataclass, field
from os import PathLike
from typing import Any, IO, BinaryIO
from zipfile import ZipFile

import msgpack

from nanover.trajectory import FrameData2
from nanover.trajectory.convert import (
    unpack_dict_frame,
    convert_dict_state_to_dictionary_change,
)

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

    def get_message_from_entry(self, entry: RecordingIndexEntry) -> dict[str, Any]:
        return msgpack.unpackb(entry.read_from(self.messagesfile))

    def iter_messages(self):
        for entry in self:
            yield entry.metadata["timestamp"], self.get_message_from_entry(entry)

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


def iter_recording_file(path: PathLike[str]):
    reader = MessageZipReader.from_path(path)
    for entry in reader:
        frame, update = None, None
        message = reader.get_message_from_entry(entry)
        if "frame" in message:
            frame = FrameData2(unpack_dict_frame(message["frame"]))
        if "state" in message:
            update = convert_dict_state_to_dictionary_change(message["state"])

        if frame is not None or update is not None:
            yield entry.metadata["timestamp"], frame, update


def parse_index(zipfile: ZipFile):
    with zipfile.open(RECORDING_INDEX_FILENAME) as index_file:
        return [
            RecordingIndexEntry(**entry) for entry in msgpack.unpackb(index_file.read())
        ]
