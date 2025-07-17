import math
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
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
    prev_frame: FrameData
    next_frame: FrameData

    @classmethod
    def make_empty(cls):
        return cls(
            timestamp=0,
            message=GetFrameResponse(),
            prev_frame=FrameData(),
            next_frame=FrameData(),
        )


@dataclass(kw_only=True)
class StateRecordingEvent:
    timestamp: int
    message: StateUpdate
    prev_state: dict
    next_state: dict

    @classmethod
    def make_empty(cls):
        return cls(timestamp=0, message=StateUpdate(), prev_state={}, next_state={})


@dataclass(kw_only=True)
class RecordingEvent:
    timestamp: int
    prev_frame_event: FrameRecordingEvent
    prev_state_event: StateRecordingEvent
    next_frame_event: Optional[FrameRecordingEvent]
    next_state_event: Optional[StateRecordingEvent]

    @property
    def prev_frame(self):
        return self.prev_frame_event.next_frame

    @property
    def next_frame(self):
        return (
            self.next_frame_event.next_frame
            if self.next_frame_event is not None
            else self.prev_frame
        )

    @property
    def prev_state(self):
        return self.prev_state_event.next_state

    @property
    def next_state(self):
        return (
            self.next_state_event.next_state
            if self.next_state_event is not None
            else self.prev_state
        )


def split_on_frame_reset(event: RecordingEvent):
    return (
        event.next_frame_event is not None
        and event.next_frame_event.message.frame_index == 0
    )


def split_on_sim_counter_change(event: RecordingEvent):
    return (
        event.next_frame_event is not None
        and SIMULATION_COUNTER in event.next_frame_event.message.frame.values
    )


def name_basic(
    *,
    input_stem: str,
    index: int,
    simulation_name: str,
):
    return f"{input_stem}-{index:02d}-{simulation_name}"


def split_recording(
    *,
    traj: Optional[PathLike[str]] = None,
    state: Optional[PathLike[str]] = None,
    split_predicate: Callable[[RecordingEvent], bool] = split_on_frame_reset,
    name_template=name_basic,
):
    """
    :param traj: Path of a NanoVer trajectory recording
    :param state: Path of a NanoVer state recording
    :param split_predicate: Predicate function that takes information about the next frame and returns True if the
    recording should split at this point.
    :param name_template: Function that generates a filename for each recording.
    """
    split_count = 0
    last_name = "UNKNOWN"
    current_base_timestamp = 0
    current_traj_out: Optional[BinaryIO] = None
    current_state_out: Optional[BinaryIO] = None

    input_path = Path(traj if traj is not None else state)
    temp_stem = f"{input_path.stem}--TEMP"

    def close_all():
        split_stem = name_template(
            input_stem=input_path.stem,
            index=split_count,
            simulation_name=last_name,
        )

        if current_traj_out is not None:
            current_traj_out.close()
            Path(f"{temp_stem}.traj").rename(input_path.parent / f"{split_stem}.traj")
        if current_state_out is not None:
            current_state_out.close()
            Path(f"{temp_stem}.state").rename(input_path.parent / f"{split_stem}.state")

    def open_all():
        nonlocal current_traj_out, current_state_out
        if traj is not None:
            current_traj_out = open(f"{temp_stem}.traj", "wb")
            write_header(current_traj_out)
        if state is not None:
            current_state_out = open(f"{temp_stem}.state", "wb")
            write_header(current_state_out)

    def dict_to_state_update(dict):
        change = DictionaryChange(updates=dict)
        state = dictionary_change_to_state_update(change)
        return state

    try:
        open_all()

        for event in iter_recording_max(traj=traj, state=state):
            last_name = event.prev_state.get("puppeteer.simulation-name", "UNKNOWN")

            if split_predicate(event):
                # end the previous recording
                close_all()

                # start the next recording
                split_count += 1
                current_base_timestamp = event.timestamp
                open_all()

                # start the next recording with the outcome of the current event
                write_entry(
                    current_traj_out,
                    0,
                    GetFrameResponse(frame_index=0, frame=event.next_frame.raw),
                )
                write_entry(
                    current_state_out, 0, dict_to_state_update(event.next_state)
                )
            else:
                timestamp = event.timestamp - current_base_timestamp
                if event.next_frame_event is not None:
                    write_entry(
                        current_traj_out, timestamp, event.next_frame_event.message
                    )
                if event.next_state_event is not None:
                    write_entry(
                        current_state_out, timestamp, event.next_state_event.message
                    )
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

    prev_frame_event = FrameRecordingEvent.make_empty()
    prev_state_event = StateRecordingEvent.make_empty()

    while next_frame is not None or next_state is not None:
        next_frame_event: Optional[FrameRecordingEvent] = None
        next_state_event: Optional[StateRecordingEvent] = None

        time = min(get_time(next_frame), get_time(next_state))
        if next_frame is not None and get_time(next_frame) == time:
            next_frame_event = next_frame
            next_frame = next(frames, None)
        if next_state is not None and get_time(next_state) == time:
            next_state_event = next_state
            next_state = next(states, None)

        yield RecordingEvent(
            timestamp=time,
            next_frame_event=next_frame_event,
            next_state_event=next_state_event,
            prev_frame_event=prev_frame_event,
            prev_state_event=prev_state_event,
        )

        prev_frame_event = next_frame_event or prev_frame_event
        prev_state_event = next_state_event or prev_state_event


def iter_frame_file_full(path: PathLike[str]):
    prev_frame = FrameData()
    with MessageRecordingReader.from_path(path) as reader:
        for entry in reader:
            message = buffer_to_frame_message(entry.buffer)
            frame_reset = message.frame_index == 0
            if frame_reset:
                prev_frame = FrameData()
            next_frame = FrameData()
            next_frame.raw.MergeFrom(prev_frame.raw)
            next_frame.raw.MergeFrom(message.frame)

            yield FrameRecordingEvent(
                timestamp=entry.timestamp,
                message=message,
                prev_frame=prev_frame,
                next_frame=next_frame,
            )

            prev_frame = next_frame


def iter_state_file_full(path: PathLike[str]):
    current_state = StateDictionary()
    prev_state = current_state.copy_content()
    with MessageRecordingReader.from_path(path) as reader:
        for entry in reader:
            message = buffer_to_state_message(entry.buffer)
            change = state_update_to_dictionary_change(message)
            current_state.update_state(None, change)
            next_state = current_state.copy_content()

            yield StateRecordingEvent(
                timestamp=entry.timestamp,
                message=message,
                prev_state=prev_state,
                next_state=next_state,
            )

            prev_state = next_state
