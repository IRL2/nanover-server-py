import random
from io import BytesIO
from itertools import zip_longest

from nanover.protocol.trajectory import GetFrameResponse
from nanover.protocol.state import StateUpdate
from nanover.recording.reading import iter_recording_entries
from nanover.recording.writing import write_header, write_entry, record_entries
from nanover.state.state_service import dictionary_change_to_state_update
from nanover.trajectory import FrameData
from nanover.utilities.change_buffers import DictionaryChange


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


def test_reads_written_frames():
    """
    Test that a written sequence of frames is read back the same.
    """
    entries = [
        (
            i * 100000 + random.randint(0, 50000),
            GetFrameResponse(frame_index=i, frame=random_frame().raw),
        )
        for i in range(1000)
    ]

    with BytesIO() as io:
        record_entries(io, entries)

        io.seek(0)

        for a, b in zip_longest(entries, iter_recording_entries(io, GetFrameResponse)):
            assert a == b


def test_reads_written_updates():
    """
    Test that a written sequence of updates is read back the same.
    """
    entries = [
        (
            i * 100000 + random.randint(0, 50000),
            dictionary_change_to_state_update(random_change()),
        )
        for i in range(1000)
    ]

    with BytesIO() as io:
        write_header(io)
        for timestamp, message in entries:
            write_entry(io, timestamp, message)

        io.seek(0)

        for a, b in zip_longest(entries, iter_recording_entries(io, StateUpdate)):
            assert a == b
