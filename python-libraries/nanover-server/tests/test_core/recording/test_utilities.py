from itertools import pairwise

import pytest

from nanover.recording.utilities import (
    iter_frame_file_full,
    iter_state_file_full,
)
from nanover.recording2.utilities import iter_recording_max
from nanover.state.state_dictionary import StateDictionary
from nanover.state.state_service import state_update_to_dictionary_change
from nanover.testing.utilities import simplify_numpy
from nanover.trajectory import FrameData2
from nanover.trajectory.convert import convert_GetFrameResponse_to_framedata2
from nanover.utilities.change_buffers import DictionaryChange

from .test_reading import (
    RECORDING_PATH,
    RECORDING_PATH_SWITCHING,
)


@pytest.mark.parametrize("path", (RECORDING_PATH, RECORDING_PATH_SWITCHING))
def test_iter_frame_full_merging(path):
    """
    Test that frame resets result in a next frame that is equal to the frame message and that otherwise the next frame
    is equal to the previous frame merged with the message frame.
    """
    for event in iter_frame_file_full(path):
        if event.message.frame_index == 0:
            frame = convert_GetFrameResponse_to_framedata2(event.message)
            assert simplify_numpy(event.next_frame.frame_dict) == simplify_numpy(
                frame.frame_dict
            )
        else:
            frame = FrameData2()
            frame.update(event.prev_frame)
            frame.update(convert_GetFrameResponse_to_framedata2(event.message))
            assert simplify_numpy(frame.frame_dict) == simplify_numpy(
                event.next_frame.frame_dict
            )


@pytest.mark.parametrize("path", (RECORDING_PATH, RECORDING_PATH_SWITCHING))
def test_iter_frame_full_sequential(path):
    """
    Test that each message's prev_frame matches the content of the previous message's next_frame.
    """
    for prev_event, next_event in pairwise(iter_frame_file_full(path)):
        assert simplify_numpy(prev_event.next_frame.frame_dict) == simplify_numpy(
            next_event.prev_frame.frame_dict
        )


@pytest.mark.parametrize("path", (RECORDING_PATH_STATE, RECORDING_PATH_STATE_SWITCHING))
def test_iter_state_full_merging(path):
    """
    Test that the next state is equal to the previous state with the update message applied.
    """
    for event in iter_state_file_full(path):
        change = state_update_to_dictionary_change(event.message)
        state = StateDictionary()
        state.update_state(None, DictionaryChange(updates=event.prev_state))
        state.update_state(None, change)
        assert event.next_state == state.copy_content()


@pytest.mark.parametrize("path", (RECORDING_PATH_STATE, RECORDING_PATH_STATE_SWITCHING))
def test_iter_state_full_sequential(path):
    """
    Test that each event's prev_frame matches the content of the previous event's next_frame.
    """
    for prev_event, next_event in pairwise(iter_state_file_full(path)):
        assert prev_event.next_state == next_event.prev_state


@pytest.mark.parametrize(
    "path",
    (
        RECORDING_PATH,
        RECORDING_PATH_SWITCHING,
    ),
)
def test_iter_max_sequential(path):
    """
    Test that each event's prev events match the content of the previous event's next events if it exists.
    """
    for prev_event, next_event in pairwise(iter_recording_max(path)):
        assert (
            prev_event.next_frame_event is None
            or prev_event.next_frame_event == next_event.prev_frame_event
        )
        assert (
            prev_event.next_state_event is None
            or prev_event.next_state_event == next_event.prev_state_event
        )


@pytest.mark.parametrize(
    "path",
    (
        RECORDING_PATH,
        RECORDING_PATH_SWITCHING,
    ),
)
def test_iter_max_contiguous(path):
    """
    Test that there is always one or both of a next frame event or next state event for each event
    """
    for event in iter_recording_max(path):
        assert event.next_frame_event is not None or event.next_state_event is not None
