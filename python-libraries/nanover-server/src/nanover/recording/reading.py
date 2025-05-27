import math
from contextlib import suppress
from dataclasses import dataclass
from itertools import groupby
from os import PathLike, SEEK_CUR
from typing import (
    Tuple,
    Iterable,
    BinaryIO,
    Protocol,
    Callable,
    TypeVar,
    Optional,
    Iterator,
)

from nanover.protocol.trajectory import GetFrameResponse
from nanover.protocol.state import StateUpdate
from nanover.state.state_dictionary import StateDictionary
from nanover.state.state_service import state_update_to_dictionary_change
from nanover.trajectory import FrameData, MissingDataError
from nanover.trajectory.frame_data import SIMULATION_COUNTER
from nanover.utilities.change_buffers import DictionaryChange

MAGIC_NUMBER = 6661355757386708963

FrameEntry = Tuple[int, int, FrameData]


@dataclass(kw_only=True)
class RecordingFileEntry:
    offset: int
    timestamp: float
    buffer: bytes


class MessageRecordingReader:
    @classmethod
    def from_path(cls, path: PathLike[str]):
        return cls(open(path, "rb"))

    def __init__(self, io: BinaryIO):
        self.io = io
        self.message_offsets: list[int] = []
        self.reindex()

    def reindex(self):
        self.io.seek(0)
        read_header(self.io)
        self.message_offsets = [entry.offset for entry in skip_buffers(self.io)]

    def close(self):
        self.io.close()

    def get_entry_at_index(self, index):
        return self.get_entry_at_offset(self.message_offsets[index])

    def get_entry_at_offset(self, offset):
        self.io.seek(offset)
        timestamp, buffer = read_buffer(self.io)
        return RecordingFileEntry(offset=offset, timestamp=timestamp, buffer=buffer)

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


def split_by_simulation_counter(
    *, traj: PathLike[str], state: Optional[PathLike[str]] = None
):
    """
    Split a trajectory recording (and optionally a corresponding state recording) into sequences that share the same
    simulation counter value.
    """

    def get_simulation_counter(triplet):
        _, frame, _ = triplet
        return frame.values.get(SIMULATION_COUNTER)

    full_view = iter_full_view(traj=traj, state=state)
    sessions = [
        list(session)
        for simulation_counter, session in groupby(
            full_view, key=get_simulation_counter
        )
        if simulation_counter is not None
    ]

    return sessions


def iter_full_view(
    *, traj: Optional[PathLike[str]] = None, state: Optional[PathLike[str]] = None
):
    """
    Iterate one or both of trajectory and state recording files, yield a timestamp and copies of both the current
    aggregate FrameData and the current aggregate state dictionary.
    """
    full_frame = FrameData()
    full_state = StateDictionary()

    for time, frame, update in iter_recording_files(traj=traj, state=state):
        if frame is not None:
            frame_reset = frame.values["index"] == 0
            if frame_reset:
                full_frame = FrameData()
            full_frame.raw.MergeFrom(frame.raw)

        if update is not None:
            full_state.update_state(None, update)

        yield time, full_frame.copy(), full_state.copy_content()


def iter_recording_files(
    *, traj: Optional[PathLike[str]] = None, state: Optional[PathLike[str]] = None
):
    """
    Iterate one or both of trajectory and state recording files, yield a timestamp and one or both of frame and update
    that occurred at that instant. Frame index is included in frame data under the key "index".
    """
    frames: Iterator[Tuple[int, int, FrameData]] = (
        iter([]) if traj is None else iter_trajectory_file(traj)
    )
    updates: Iterator[Tuple[int, StateUpdate]] = (
        iter([]) if state is None else iter_state_file(state)
    )

    next_frame = next(frames, None)
    next_update = next(updates, None)

    def get_time(entry):
        return math.inf if entry is None else entry[0]

    while next_frame is not None or next_update is not None:
        frame: Optional[FrameData] = None
        update: Optional[DictionaryChange] = None

        time = min(get_time(next_frame), get_time(next_update))
        if next_frame is not None and get_time(next_frame) == time:
            _, index, frame = next_frame
            frame.values["index"] = index
            next_frame = next(frames, None)
        if next_update is not None and get_time(next_update) == time:
            _, update = next_update
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
        yield elapsed, buffer_to_message(buffer, message_type)


def iter_recording_buffers(io: BinaryIO):
    """
    Iterate over elements of a recording, yielding pairs of a timestamp in microseconds and a buffer of bytes.

    :param io: Stream of bytes to read.
    """
    read_header(io)
    yield from ((entry.timestamp, entry.buffer) for entry in iter_buffers(io))


def iter_trajectory_recording(io: BinaryIO):
    for elapsed, get_frame_response in iter_recording_entries(io, GetFrameResponse):
        frame_index: int = get_frame_response.frame_index
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


def buffer_to_frame_message(buffer):
    return buffer_to_message(buffer, GetFrameResponse)


def buffer_to_state_message(buffer):
    return buffer_to_message(buffer, StateUpdate)


def buffer_to_message(buffer, message_type: Callable[[], TMessage]):
    instance = message_type()
    instance.ParseFromString(buffer)
    return instance


def iter_buffers(io: BinaryIO):
    with suppress(EOFError):
        while True:
            position = io.tell()
            timestamp, buffer = read_buffer(io)
            yield RecordingFileEntry(
                offset=position, timestamp=timestamp, buffer=buffer
            )


def skip_buffers(io: BinaryIO):
    with suppress(EOFError):
        while True:
            position = io.tell()
            timestamp, _ = skip_buffer(io)
            yield RecordingFileEntry(
                offset=position, timestamp=timestamp, buffer=bytes()
            )


def read_header(io: BinaryIO):
    supported_format_versions = (2,)

    magic_number = read_u64(io)
    if magic_number != MAGIC_NUMBER:
        raise InvalidMagicNumber

    format_version = read_u64(io)
    if format_version not in supported_format_versions:
        raise UnsupportedFormatVersion(format_version, supported_format_versions)


def read_buffer(io: BinaryIO):
    timestamp = read_u128(io)
    record_size = read_u64(io)
    buffer = io.read(record_size)
    return timestamp, buffer


def skip_buffer(io: BinaryIO):
    timestamp = read_u128(io)
    record_size = read_u64(io)
    io.seek(record_size, SEEK_CUR)
    return timestamp, record_size


def read_u64(io: BinaryIO):
    return int.from_bytes(read(io, 8), "little", signed=False)


def read_u128(io: BinaryIO):
    return int.from_bytes(read(io, 16), "little", signed=False)


def read(io: BinaryIO, count: int):
    buffer = io.read(count)
    if len(buffer) < count:
        raise EOFError
    return buffer
