from typing import Any

from nanover.recording.old.utilities import iter_recording_max
from nanover.recording.reading import MessageEvent
from nanover.recording.writing import record_messages
from nanover.state.state_service import state_update_to_dictionary_change
from .frame_data import FrameData as FrameDataOld
from nanover.trajectory.keys import FRAME_INDEX
from nanover.trajectory.convert import converters, pack_identity, pack_force_list


def convert_old_recording(out_path, *, traj, state):
    record_messages(out_path, message_events_from_old_recording(traj=traj, state=state))


def message_events_from_old_recording(*, traj, state):
    for event in iter_recording_max(traj=traj, state=state):
        if event.next_frame_event is not None:
            frame_response = event.next_frame_event.message
            frame = pack_grpc_frame(FrameDataOld(frame_response.frame))
            frame[FRAME_INDEX] = frame_response.frame_index
            message = {"frame": frame}
            yield MessageEvent(
                timestamp=event.timestamp,
                message=message,
            )
        if event.next_state_event is not None:
            change = state_update_to_dictionary_change(event.next_state_event.message)
            message = {
                "state": {"updates": change.updates, "removals": list(change.removals)}
            }
            yield MessageEvent(
                timestamp=event.timestamp,
                message=message,
            )


def pack_grpc_frame(frame: FrameDataOld) -> dict[str, Any]:
    data = {}

    for key in frame.value_keys:
        converter = converters.get(key, pack_identity)
        data[key] = converter.pack(frame.values[key])

    for key in frame.array_keys:
        converter = converters.get(key, pack_force_list)
        data[key] = converter.pack(frame.arrays[key])

    return data
