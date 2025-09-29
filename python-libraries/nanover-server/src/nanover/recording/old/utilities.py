import math
from dataclasses import dataclass
from os import PathLike
from typing import Any

from .protocol.trajectory import GetFrameResponse
from .protocol.state import StateUpdate
from nanover.recording.old.reading import (
    MessageRecordingReader,
    buffer_to_frame_message,
    buffer_to_state_message,
)
from nanover.state.state_dictionary import StateDictionary
from nanover.trajectory import FrameData
from .frame_data import FrameData as FrameDataOld
from ...trajectory.convert import (
    unpack_dict_frame,
    pack_identity,
    pack_force_list,
    converters,
)
from ...utilities.change_buffers import DictionaryChange
from ...utilities.protobuf_utilities import struct_to_dict


def convert_grpc_frame_to_dict_frame(grpc_frame) -> dict[str, Any]:
    return unpack_dict_frame(pack_grpc_frame(grpc_frame))


def convert_GetFrameResponse_to_framedata(response: GetFrameResponse) -> FrameData:
    frame = FrameData(convert_grpc_frame_to_dict_frame(FrameDataOld(response.frame)))
    frame.frame_index = response.frame_index
    return frame


def pack_grpc_frame(frame: FrameDataOld) -> dict[str, Any]:
    data = {}

    for key in frame.value_keys:
        converter = converters.get(key, pack_identity)
        data[key] = converter.pack(frame.values[key])

    for key in frame.array_keys:
        converter = converters.get(key, pack_force_list)
        data[key] = converter.pack(frame.arrays[key])

    return data


def state_update_to_dictionary_change(update: StateUpdate) -> DictionaryChange:
    """
    Convert a protobuf StateUpdate to a DictionaryChange.

    :param update: a protobuf StateUpdate which encodes key removals as keys
        with a protobuf null value.
    :return: an equivalent DictionaryChange representing the key changes and
        key removals of the StateUpdate.
    """
    removals = set()
    updates = {}

    for key, value in struct_to_dict(update.changed_keys).items():
        if value is None:
            removals.add(key)
        else:
            updates[key] = value

    return DictionaryChange(removals=removals, updates=updates)


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
    next_frame_event: FrameRecordingEvent | None
    next_state_event: StateRecordingEvent | None

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


def iter_recording_max(
    *, traj: PathLike[str] | None = None, state: PathLike[str] | None = None
):
    """
    Iterate recording file(s) yielding recording events in timestamp order, with each event containing the full
    information of previous frame, previous state, current message, next frame, and next state.
    """
    frames = iter([]) if traj is None else iter_frame_file_full(traj)
    states = iter([]) if state is None else iter_state_file_full(state)

    next_frame = next(frames, None)
    next_state = next(states, None)

    def get_time(entry):
        return math.inf if entry is None else entry.timestamp

    prev_frame_event = FrameRecordingEvent.make_empty()
    prev_state_event = StateRecordingEvent.make_empty()

    while next_frame is not None or next_state is not None:
        next_frame_event: FrameRecordingEvent | None = None
        next_state_event: StateRecordingEvent | None = None

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
    """
    Yield an event for every message in the recording that contains the timestamp, message, and frame data before and
    after the message was applied.
    """
    prev_frame = FrameData()
    with MessageRecordingReader.from_path(path) as reader:
        for entry in reader:
            message = buffer_to_frame_message(entry.buffer)
            frame = convert_GetFrameResponse_to_framedata(message)

            next_frame = prev_frame.copy()
            next_frame.update(frame)

            yield FrameRecordingEvent(
                timestamp=entry.timestamp,
                message=message,
                prev_frame=prev_frame,
                next_frame=next_frame,
            )

            prev_frame = next_frame


def iter_state_file_full(path: PathLike[str]):
    """
    Yield an event for every message in the recording that contains the timestamp, message, and state before and after
    the message was applied.
    """
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
