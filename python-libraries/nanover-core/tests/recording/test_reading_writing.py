import random
from io import BytesIO
from itertools import zip_longest

import pytest

from nanover.protocol.trajectory import GetFrameResponse
from nanover.protocol.state import StateUpdate
from nanover.recording.reading import iter_recording_entries
from nanover.recording.writing import write_header, write_entry, record_entries
from nanover.state.state_service import dictionary_change_to_state_update
from nanover.trajectory import FrameData
from nanover.utilities.change_buffers import DictionaryChange


def random_frame_message():
    return GetFrameResponse(
        frame_index=random.randint(0, 1000), frame=random_frame().raw
    )


def random_state_message():
    return dictionary_change_to_state_update(random_change())


MESSAGE_TYPE_FACTORY_PAIRS = (
    (GetFrameResponse, random_frame_message),
    (StateUpdate, random_state_message),
)


def random_frame():
    frame = FrameData()
    frame.particle_count = random.randint(10, 100)
    frame.particle_positions = [
        [random.random(), random.random(), random.random()]
        for _ in range(frame.particle_count)
    ]
    return frame


def random_change():
    return DictionaryChange(
        updates={
            "a" * random.randint(3, 8): random.random()
            for _ in range(random.randint(5, 10))
        },
        removals={
            "a" * random.randint(3, 8): random.random()
            for _ in range(random.randint(5, 10))
        },
    )


@pytest.mark.parametrize("pair", MESSAGE_TYPE_FACTORY_PAIRS)
def test_reads_written_messages(pair):
    """
    Test that a written sequence of messages is read back the same.
    """
    message_type, random_message = pair
    entries = [
        (
            i * 100000 + random.randint(0, 50000),
            random_message(),
        )
        for i in range(1000)
    ]

    with BytesIO() as io:
        record_entries(io, entries)

        io.seek(0)

        for a, b in zip_longest(entries, iter_recording_entries(io, message_type)):
            assert a == b
