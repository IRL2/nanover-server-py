from time import perf_counter_ns
from typing import BinaryIO, Protocol, Optional, Iterable

from nanover.recording.reading import MAGIC_NUMBER

WRITE_VERSION = 2


class Serializable(Protocol):
    def SerializeToString(self) -> bytes: ...


def record_messages_to_file(path, messages: Iterable[Serializable]):
    """
    Write a sequence of messages to a recording file.

    :param path: Path of recording file to write to.
    :param messages: Iterable sequence of messages to record.
    """
    with open(path, "wb") as outfile:
        yield from record_messages(outfile, messages)


def record_messages(io: BinaryIO, messages, start_time: Optional[int] = None):
    if start_time is None:
        start_time = perf_counter_ns()

    def timestamp():
        return int((perf_counter_ns() - start_time) / 1000)

    entries = (timestamp(), message for message in messages)
    record_entries(io, entries)


def record_entries(io: BinaryIO, entries):
    write_header(io)
    for (timestamp, message) in entries:
        write_entry(io, timestamp, message)


def write_header(io: BinaryIO):
    write_u64(io, MAGIC_NUMBER)
    write_u64(io, WRITE_VERSION)


def write_entry(io: BinaryIO, timestamp: int, message: Serializable):
    buffer = message.SerializeToString()
    write_u128(io, timestamp)
    write_u64(io, len(buffer))
    io.write(buffer)


def write_u64(io: BinaryIO, value: int):
    io.write(value.to_bytes(8, "little", signed=False))


def write_u128(io: BinaryIO, value: int):
    io.write(value.to_bytes(16, "little", signed=False))
