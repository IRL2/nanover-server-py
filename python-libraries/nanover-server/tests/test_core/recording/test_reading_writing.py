from io import BytesIO
from itertools import zip_longest

from hypothesis import strategies as st
from hypothesis import given

from nanover.testing.strategies import packed_frame_dicts, state_updates
from nanover.recording.reading import NanoverRecordingReader
from nanover.recording.writing import record_messages, MessageEvent


@st.composite
def message_event_frames(draw):
    timestamp = draw(st.integers(min_value=0, max_value=2**32 - 1))
    message = {"frame": draw(packed_frame_dicts())}
    return MessageEvent(timestamp=timestamp, message=message)


@st.composite
def message_event_states(draw):
    timestamp = draw(st.integers(min_value=0, max_value=2**32 - 1))
    message = {"state": draw(state_updates())}
    return MessageEvent(timestamp=timestamp, message=message)


message_events_all = st.lists(st.one_of(message_event_states(), message_event_frames()))


@given(message_events=message_events_all)
def test_reads_written_messages(message_events):
    """
    Test that a written sequence of messages is read back the same.
    """
    with BytesIO() as io:
        record_messages(io, message_events)
        reader = NanoverRecordingReader.from_io(io)
        for a, b in zip_longest(message_events, reader.iter_message_events()):
            assert a == b
