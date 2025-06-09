import math
from dataclasses import dataclass
from os import PathLike
from typing import Optional, BinaryIO, Callable

from nanover.protocol.trajectory import GetFrameResponse
from nanover.protocol.state import StateUpdate
from nanover.recording.reading import (
    MessageRecordingReader,
    buffer_to_frame_message,
    buffer_to_state_message,
)
from nanover.recording.writing import write_entry, write_header
from nanover.state.state_dictionary import StateDictionary
from nanover.state.state_service import (
    dictionary_change_to_state_update,
    state_update_to_dictionary_change,
)
from nanover.trajectory import FrameData
from nanover.trajectory.frame_data import SIMULATION_COUNTER
from nanover.utilities.change_buffers import DictionaryChange


@dataclass(kw_only=True)
class FrameRecordingEvent:
    timestamp: int
    message: GetFrameResponse
    current_frame: FrameData


@dataclass(kw_only=True)
class StateRecordingEvent:
    timestamp: int
    message: StateUpdate
    current_state: dict


@dataclass(kw_only=True)
class RecordingEvent:
    timestamp: int
    frame_event: Optional[FrameRecordingEvent]
    state_event: Optional[StateRecordingEvent]


def split_on_frame_reset(event: RecordingEvent):
    return event.frame_event is not None and event.frame_event.message.frame_index == 0


def split_on_sim_counter_change(event: RecordingEvent):
    return (
        event.frame_event is not None
        and SIMULATION_COUNTER in event.frame_event.message.frame.values
    )


def split_recording(
    *,
    traj: Optional[PathLike[str]] = None,
    state: Optional[PathLike[str]] = None,
    split_predicate: Callable[[RecordingEvent], bool] = split_on_frame_reset,
):
    """
    :param traj: Path of a NanoVer trajectory recording
    :param state: Path of a NanoVer state recording
    :param split_predicate: Predicate function that takes information about the next frame and returns True if the
    recording should split at this point.
    """
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

        prev_frame = None
        prev_state = None

        for event in iter_recording_max(traj=traj, state=state):
            if event.frame_event is not None:
                prev_frame = event.frame_event
            if event.frame_event is not None:
                prev_state = event.state_event

            if split_predicate(event):
                close_all()
                split_count += 1
                current_base_timestamp = event.timestamp
                open_all()
                if prev_frame is not None:
                    frame = GetFrameResponse(
                        frame_index=0, frame=prev_frame.current_frame.raw
                    )
                    write_entry(current_traj_out, 0, frame)
                if prev_state is not None:
                    state_update = dict_to_state_update(prev_state.current_state)
                    write_entry(current_state_out, 0, state_update)
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
    with MessageRecordingReader.from_path(path) as reader:
        for entry in reader:
            message = buffer_to_frame_message(entry.buffer)
            frame_reset = message.frame_index == 0
            if frame_reset:
                current_frame = FrameData()
            current_frame.raw.MergeFrom(message.frame)

            yield FrameRecordingEvent(
                timestamp=entry.timestamp,
                message=message,
                current_frame=current_frame,
            )


def iter_state_file_full(path: PathLike[str]):
    current_state = StateDictionary()
    with MessageRecordingReader.from_path(path) as reader:
        for entry in reader:
            message = buffer_to_state_message(entry.buffer)
            change = state_update_to_dictionary_change(message)
            current_state.update_state(None, change)

            yield StateRecordingEvent(
                timestamp=entry.timestamp,
                message=message,
                current_state=current_state.copy_content(),
            )
