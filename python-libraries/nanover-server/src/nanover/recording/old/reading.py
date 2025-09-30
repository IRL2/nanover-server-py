from contextlib import suppress
from dataclasses import dataclass
from os import PathLike, SEEK_CUR
from typing import (
    BinaryIO,
    Protocol,
    Callable,
    TypeVar,
)

from nanover.protocol.trajectory import GetFrameResponse
from nanover.protocol.state import StateUpdate

MAGIC_NUMBER = 6661355757386708963


@dataclass(kw_only=True)
class RecordingFileEntry:
    offset: int
    timestamp: float
    buffer: bytes


class Parseable(Protocol):
    def ParseFromString(self, _: bytes) -> bytes: ...


class InvalidMagicNumber(Exception):
    """
    The magic number read from a file is not the expected one.

    The file may not be in the correct format, or it may be corrupted.
    """


class UnsupportedFormatVersion(Exception):
    """
    The version of the file format is not supported by the parser.
    """

    def __init__(self, format_version: int, supported_format_versions: tuple[int]):
        self._format_version = format_version
        self._supported_format_versions = supported_format_versions

    def __str__(self) -> str:
        return (
            f"Version {self._format_version} of the format is not supported "
            f"by this parser. The supported versions are "
            f"{self._supported_format_versions}."
        )


TMessage = TypeVar("TMessage", bound=Parseable)


class MessageRecordingReader:
    """
    Wraps a NanoVer recording of gRPC messages to provide fast and convenient random access.
    """

    @classmethod
    def from_path(cls, path: PathLike[str]):
        """
        Read a recording from a filepath.
        """
        return cls.from_io(open(path, "rb"))

    @classmethod
    def from_io(cls, io: BinaryIO):
        """
        Read a recording from binary data source.
        """
        reader = cls(io)
        reader.reindex()
        return reader

    def __init__(self, io: BinaryIO):
        self.io = io
        self.message_offsets: list[int] = []

    def reindex(self):
        """
        Read the file from beginning to end and recording the byte position of each message for later access.
        """
        self.io.seek(0)
        read_header(self.io)
        self.message_offsets = [entry.offset for entry in skip_buffers(self.io)]

    def close(self):
        self.io.close()

    def get_entry_at_index(self, index):
        """
        Return an entry for the nth message in the file.

        :param index: Index of the message in the sequence of all messages in the data.
        """
        return self.get_entry_at_offset(self.message_offsets[index])

    def get_entry_at_offset(self, offset):
        """
        Read and return a message entry from a specific byte offset in the data.

        :param offset: Offset in bytes to begin reading a message entry.
        """
        self.io.seek(offset)
        timestamp, buffer = read_buffer(self.io)
        return RecordingFileEntry(offset=offset, timestamp=timestamp, buffer=buffer)

    def iter_messages(self, message_type: Callable[[], TMessage]):
        for entry in self:
            yield (entry.timestamp, buffer_to_message(entry.buffer, message_type))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        for offset in self.message_offsets:
            yield self.get_entry_at_offset(offset)

    def __len__(self):
        return len(self.message_offsets)

    def __getitem__(self, index):
        return self.get_entry_at_index(index)


def buffer_to_frame_message(buffer):
    return buffer_to_message(buffer, GetFrameResponse)


def buffer_to_state_message(buffer):
    return buffer_to_message(buffer, StateUpdate)


def buffer_to_message(buffer, message_type: Callable[[], TMessage]):
    instance = message_type()
    instance.ParseFromString(buffer)
    return instance


def skip_buffers(io: BinaryIO):
    with suppress(EOFError):
        while True:
            position = io.tell()
            timestamp, _ = skip_buffer(io)
            yield RecordingFileEntry(
                offset=position, timestamp=timestamp, buffer=bytes()
            )


def read_buffer(io: BinaryIO):
    timestamp, buffer_size = read_entry_info(io)
    buffer = io.read(buffer_size)
    return timestamp, buffer


def skip_buffer(io: BinaryIO):
    timestamp, buffer_size = read_entry_info(io)
    io.seek(buffer_size, SEEK_CUR)
    return timestamp, buffer_size


def read_entry_info(io: BinaryIO):
    timestamp = read_u128(io)
    buffer_size = read_u64(io)
    return timestamp, buffer_size


def read_header(io: BinaryIO):
    supported_format_versions = (2,)

    magic_number = read_u64(io)
    if magic_number != MAGIC_NUMBER:
        raise InvalidMagicNumber

    format_version = read_u64(io)
    if format_version not in supported_format_versions:
        raise UnsupportedFormatVersion(format_version, supported_format_versions)

    return format_version


def read_u64(io: BinaryIO):
    return int.from_bytes(read(io, 8), "little", signed=False)


def read_u128(io: BinaryIO):
    return int.from_bytes(read(io, 16), "little", signed=False)


def read(io: BinaryIO, count: int):
    buffer = io.read(count)
    if len(buffer) < count:
        raise EOFError
    return buffer
