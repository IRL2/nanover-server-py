from os import PathLike

from nanover.recording.utilities import (
    RecordingEvent,
    FrameRecordingEvent,
    StateRecordingEvent,
)
from nanover.recording2.reading import iter_recording_file


def iter_recording_max(path: PathLike[str]):
    """
    Iterate recording file yielding recording events in timestamp order, with each event containing the full
    information of previous frame, previous state, current message, next frame, and next state.
    """
    entries = iter_recording_file(path)

    next_entry = next(entries, None)

    prev_frame_event = FrameRecordingEvent.make_empty()
    prev_state_event = StateRecordingEvent.make_empty()

    for timestamp, frame, update in iter_recording_file(path):
        next_frame_event = (
            FrameRecordingEvent(
                timestamp=timestamp,
            )
            if frame is not None
            else None
        )

        yield RecordingEvent(
            timestamp=timestamp,
        )
