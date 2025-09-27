from nanover.recording.utilities import iter_recording_max
from nanover.recording2.reading import MessageEvent
from nanover.recording2.writing import record_messages
from nanover.state.state_service import state_update_to_dictionary_change
from nanover.trajectory import FrameData
from nanover.trajectory.frame_data import FRAME_INDEX
from nanover.trajectory.convert import pack_grpc_frame


def convert_old_recording(out_path, *, traj, state):
    record_messages(out_path, message_events_from_old_recording(traj=traj, state=state))


def message_events_from_old_recording(*, traj, state):
    for event in iter_recording_max(traj=traj, state=state):
        if event.next_frame_event is not None:
            frame_response = event.next_frame_event.message
            frame = pack_grpc_frame(FrameData(frame_response.frame))
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
