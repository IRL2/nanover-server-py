import math
from dataclasses import dataclass, KW_ONLY
from os import PathLike
from typing import Optional, BinaryIO, Callable

from nanover.protocol.trajectory import GetFrameResponse
from nanover.protocol.state import StateUpdate

from nanover.recording.reading import iter_recording_entries
from nanover.recording.writing import write_header, write_entry
from nanover.state.state_dictionary import StateDictionary
from nanover.state.state_service import (
    state_update_to_dictionary_change,
    dictionary_change_to_state_update,
)
from nanover.trajectory import FrameData
from nanover.utilities.change_buffers import DictionaryChange


@dataclass
class FrameRecordingEvent:
    _: KW_ONLY
    timestamp: int
    message: GetFrameResponse
    current_frame: FrameData


@dataclass
class StateRecordingEvent:
    _: KW_ONLY
    timestamp: int
    message: StateUpdate
    current_state: dict


@dataclass
class RecordingEvent:
    _: KW_ONLY
    timestamp: int
    frame_event: Optional[FrameRecordingEvent]
    state_event: Optional[StateRecordingEvent]


def split_recording(
    *,
    traj: Optional[PathLike[str]] = None,
    state: Optional[PathLike[str]] = None,
    split_predicate=Callable[[RecordingEvent], bool],
):
    split_count = 0
    current_base_timestamp = 0
    current_traj_out: Optional[BinaryIO] = None
    current_state_out: Optional[BinaryIO] = None

    def close_all():
        if current_state_out is not None:
            current_traj_out.close()
        if current_state_out is not None:
            current_state_out.close()

    def open_all():
        nonlocal current_traj_out, current_state_out
        if traj is not None:
            current_traj_out = open(f"{traj}--SPLIT--{split_count}.traj", "wb")
            write_header(current_traj_out)
        if state is not None:
            current_state_out = open(f"{state}--SPLIT--{split_count}.state", "wb")
            write_header(current_state_out)

    def dict_to_state_update(dict):
        change = DictionaryChange(updates=dict)
        state = dictionary_change_to_state_update(change)
        return state

    try:
        open_all()

        for event in iter_recording_max(traj=traj, state=state):
            if split_predicate(event):
                close_all()
                split_count += 1
                current_base_timestamp = event.timestamp
                open_all()
                if event.frame_event is not None:
                    frame = GetFrameResponse(
                        frame_index=0, frame=event.frame_event.current_frame.raw
                    )
                    write_entry(current_traj_out, 0, frame)
                if event.state_event is not None:
                    state = dict_to_state_update(event.state_event.current_state)
                    write_entry(current_state_out, 0, state)
            else:
                timestamp = event.timestamp - current_base_timestamp
                if event.frame_event is not None:
                    write_entry(current_traj_out, timestamp, event.frame_event.message)
                if event.state_event is not None:
                    write_entry(current_state_out, timestamp, event.state_event.message)

    finally:
        close_all()


def iter_recording_max(
    *, traj: Optional[PathLike[str]] = None, state: Optional[PathLike[str]] = None
):
    frames = iter([]) if traj is None else iter_frame_file_full(traj)
    states = iter([]) if state is None else iter_state_file_full(state)

    next_frame = next(frames, None)
    next_state = next(states, None)

    def get_time(entry):
        return math.inf if entry is None else entry.timestamp

    while next_frame is not None or next_state is not None:
        frame: Optional[FrameRecordingEvent] = None
        state: Optional[StateRecordingEvent] = None

        time = min(get_time(next_frame), get_time(next_state))
        if next_frame is not None and get_time(next_frame) == time:
            frame = next_frame
            next_frame = next(frames, None)
        if next_state is not None and get_time(next_state) == time:
            state = next_state
            next_state = next(states, None)

        yield RecordingEvent(
            timestamp=time,
            frame_event=frame,
            state_event=state,
        )


def iter_frame_file_full(path: PathLike[str]):
    current_frame = FrameData()

    for timestamp, message in iter_from_file(path, iter_trajectory_recording):
        frame_reset = message.frame_index == 0
        if frame_reset:
            current_frame = FrameData()
        current_frame.raw.MergeFrom(message.frame)

        yield FrameRecordingEvent(
            timestamp=timestamp,
            message=message,
            current_frame=current_frame,
        )


def iter_state_file_full(path: PathLike[str]):
    current_state = StateDictionary()

    for timestamp, message in iter_from_file(path, iter_state_recording):
        change = state_update_to_dictionary_change(message)
        current_state.update_state(None, change)

        yield StateRecordingEvent(
            timestamp=timestamp,
            message=message,
            current_state=current_state.copy_content(),
        )


def iter_trajectory_recording(io: BinaryIO):
    yield from iter_recording_entries(io, GetFrameResponse)


def iter_state_recording(io: BinaryIO):
    yield from iter_recording_entries(io, StateUpdate)


def iter_from_file(path, func):
    with open(path, "rb") as infile:
        yield from func(infile)
