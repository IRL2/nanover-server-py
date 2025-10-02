from itertools import pairwise
import pytest
from nanover.recording.reading import NanoverRecordingReader
from .test_reading import RECORDING_PATH, RECORDING_PATH_SWITCHING


@pytest.mark.parametrize("path", (RECORDING_PATH, RECORDING_PATH_SWITCHING))
def test_iter_max_sequential(path):
    """
    Test that each event's prev events match the content of the previous event's next events if it exists.
    """
    with NanoverRecordingReader.from_path(path) as reader:
        for prev_event, next_event in pairwise(reader.iter_max()):
            assert (
                prev_event.next_frame_event is None
                or prev_event.next_frame_event == next_event.prev_frame_event
            )
            assert (
                prev_event.next_state_event is None
                or prev_event.next_state_event == next_event.prev_state_event
            )


@pytest.mark.parametrize("path", (RECORDING_PATH, RECORDING_PATH_SWITCHING))
def test_iter_max_contiguous(path):
    """
    Test that there is always one or both of a next frame event or next state event for each event
    """
    with NanoverRecordingReader.from_path(path) as reader:
        for event in reader.iter_max():
            assert (
                event.next_frame_event is not None or event.next_state_event is not None
            )
