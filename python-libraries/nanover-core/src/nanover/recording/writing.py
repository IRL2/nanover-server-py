from typing import BinaryIO, Protocol

from nanover.recording.reading import MAGIC_NUMBER

WRITE_VERSION = 2


class Serializable(Protocol):
    def SerializeToString(self) -> bytes: ...


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
