from typing import Tuple, Iterable, BinaryIO, Protocol, Callable, TypeVar, Optional
from nanover.protocol.trajectory import GetFrameResponse
from nanover.protocol.state import StateUpdate
from nanover.state.state_service import state_update_to_dictionary_change
from nanover.trajectory import FrameData, MissingDataError

MAGIC_NUMBER = 6661355757386708963

FrameEntry = Tuple[int, int, FrameData]


def iter_recording_files(*, traj: Optional[str] = None, state: Optional[str] = None):
    """
    Iterate one or both of trajectory and state recording files, yield a timestamp and one or both of frame and update
    that occurred at that instant.
    """
    frames = iter([]) if traj is None else iter_trajectory_file(traj)
    updates = iter([]) if state is None else iter_state_file(state)

    next_frame = next(frames, None)
    next_update = next(updates, None)

    while next_frame is not None and next_update is not None:
        frame = None
        update = None

        time = min(next_frame[0], next_update[0])
        if next_frame is not None and next_frame[0] == time:
            frame = next_frame[2]
            next_frame = next(frames, None)
        if next_update is not None and next_update[0] == time:
            update = next_update[1]
            next_update = next(updates, None)

        yield time, frame, update


def iter_trajectory_file(path):
    """
    Iterate over all frame updates in a recording file.

    :param path: Path of recording file to read from.
    """
    with open(path, "rb") as infile:
        yield from iter_trajectory_recording(infile)


def iter_state_file(path):
    """
    Iterate over all state updates in a recording file.

    :param path: Path of recording file to read from.
    """
    with open(path, "rb") as infile:
        yield from iter_state_recording(infile)


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


class Parseable(Protocol):
    def ParseFromString(self, _: bytes) -> bytes: ...


TMessage = TypeVar("TMessage", bound=Parseable)


def iter_recording_entries(io: BinaryIO, message_type: Callable[[], TMessage]):
    for elapsed, buffer in iter_recording_buffers(io):
        instance = message_type()
        instance.ParseFromString(buffer)
        yield elapsed, instance


def iter_recording_buffers(io: BinaryIO):
    """
    Iterate over elements of a recording, yielding pairs of a timestamp in microseconds and a buffer of bytes.

    :param io: Stream of bytes to read.
    """
    supported_format_versions = (2,)

    magic_number = read_u64(io)
    if magic_number != MAGIC_NUMBER:
        raise InvalidMagicNumber
    format_version = read_u64(io)
    if format_version not in supported_format_versions:
        raise UnsupportedFormatVersion(format_version, supported_format_versions)
    while True:
        try:
            elapsed = read_u128(io)
            record_size = read_u64(io)
            buffer = io.read(record_size)
        except EOFError:
            break
        yield elapsed, buffer


def iter_trajectory_recording(io: BinaryIO):
    for elapsed, get_frame_response in iter_recording_entries(io, GetFrameResponse):
        frame_index = get_frame_response.frame_index
        frame = FrameData(get_frame_response.frame)
        yield (elapsed, frame_index, frame)


def iter_state_recording(io: BinaryIO):
    for elapsed, state_update in iter_recording_entries(io, StateUpdate):
        dictionary_change = state_update_to_dictionary_change(state_update)
        yield (elapsed, dictionary_change)


def iter_trajectory_with_elapsed_integrated(frames: Iterable[FrameEntry]):
    for elapsed, frame_index, frame in frames:
        frame.values["elapsed"] = elapsed
        yield (elapsed, frame_index, frame)


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


def read_u64(io: BinaryIO):
    return int.from_bytes(read(io, 8), "little", signed=False)


def read_u128(io: BinaryIO):
    return int.from_bytes(read(io, 16), "little", signed=False)


def read(io: BinaryIO, count: int):
    buffer = io.read(count)
    if len(buffer) < count:
        raise EOFError
    return buffer
