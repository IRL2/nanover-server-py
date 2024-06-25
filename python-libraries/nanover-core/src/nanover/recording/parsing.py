from typing import Tuple, Iterable
from nanover.protocol.trajectory import GetFrameResponse
from nanover.protocol.state import StateUpdate
from nanover.state.state_service import state_update_to_dictionary_change
from nanover.trajectory import FrameData, MissingDataError

MAGIC_NUMBER = 6661355757386708963

FrameEntry = Tuple[int, int, FrameData]


class Unpacker:
    """
    Unpack data representations from a buffer of bytes.

    The unpacking methods return the requested value from the next bytes of the
    buffer and move the cursor forward. They raise an `IndexError` if the
    buffer is too short to fulfill the request.
    """

    _buffer: bytes
    _cursor: int

    def __init__(self, data: bytes):
        self._buffer = data
        self._cursor = 0

    def unpack_bytes(self, n_bytes: int) -> bytes:
        """
        Get the next `n_bytes` bytes from the buffer.

        The method raises a `ValueError` if the requested number of bytes is
        negative.
        """
        if n_bytes < 0:
            raise ValueError("Cannot unpack a negative number of bytes.")
        end = self._cursor + n_bytes
        if end >= len(self._buffer):
            raise IndexError("Not enough bytes left in the buffer.")
        bytes_to_return = self._buffer[self._cursor : end]
        self._cursor = end
        return bytes_to_return

    def unpack_u64(self) -> int:
        """
        Get an unsigned 64 bits integer from the next bytes of the buffer.
        """
        buffer = self.unpack_bytes(8)
        return int.from_bytes(buffer, "little", signed=False)

    def unpack_u128(self) -> int:
        """
        Get an unsigned 128 bits integer from the next bytes of the buffer.
        """
        buffer = self.unpack_bytes(16)
        return int.from_bytes(buffer, "little", signed=False)


class InvalidMagicNumber(Exception):
    """
    The magic number read from a file is not the expected one.

    The file may not be in the correct format, or it may be corrupted.
    """


class UnsupportedFormatVersion(Exception):
    """
    The version of the file format is not supported by hthe parser.
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


def iter_recording_buffers(unpacker: Unpacker):
    supported_format_versions = (2,)
    magic_number = unpacker.unpack_u64()
    if magic_number != MAGIC_NUMBER:
        raise InvalidMagicNumber
    format_version = unpacker.unpack_u64()
    if format_version not in supported_format_versions:
        raise UnsupportedFormatVersion(format_version, supported_format_versions)
    while True:
        try:
            elapsed = unpacker.unpack_u128()
            record_size = unpacker.unpack_u64()
            buffer = unpacker.unpack_bytes(record_size)
        except IndexError:
            break
        yield elapsed, buffer


def iter_trajectory_recording(unpacker: Unpacker):
    for elapsed, buffer in iter_recording_buffers(unpacker):
        get_frame_response = GetFrameResponse()
        get_frame_response.ParseFromString(buffer)
        frame_index = get_frame_response.frame_index
        frame = FrameData(get_frame_response.frame)
        yield (elapsed, frame_index, frame)


def iter_trajectory_with_elapsed_integrated(frames: Iterable[FrameEntry]):
    for elapsed, frame_index, frame in frames:
        frame.values["elapsed"] = elapsed
        yield (elapsed, frame_index, frame)


def iter_trajectory_file(path):
    with open(path, "rb") as infile:
        data = infile.read()
    unpacker = Unpacker(data)
    yield from iter_trajectory_recording(unpacker)


def advance_to_first_particle_frame(frames: Iterable[FrameEntry]):
    for elapsed, frame_index, frame in frames:
        try:
            particle_count = frame.particle_count
        except MissingDataError:
            pass
        else:
            if particle_count > 0:
                break
    else:
        return

    yield (elapsed, frame_index, frame)
    yield from frames


def advance_to_first_coordinate_frame(frames: Iterable[FrameEntry]):
    for elapsed, frame_index, frame in frames:
        try:
            frame.particle_positions
        except MissingDataError:
            pass
        else:
            break
    else:
        return

    yield (elapsed, frame_index, frame)
    yield from frames


def iter_state_recording(unpacker: Unpacker):
    for elapsed, buffer in iter_recording_buffers(unpacker):
        state_update = StateUpdate()
        state_update.ParseFromString(buffer)
        dictionary_change = state_update_to_dictionary_change(state_update)
        yield (elapsed, dictionary_change)


def iter_state_file(path):
    with open(path, "rb") as infile:
        data = infile.read()
    unpacker = Unpacker(data)
    yield from iter_state_recording(unpacker)
