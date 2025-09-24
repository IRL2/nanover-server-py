from pathlib import Path

import pytest

from nanover.recording.reading import (
    iter_recording_files,
    iter_full_view,
    split_by_simulation_counter,
    MessageRecordingReader,
)
from nanover.recording2.reading import MessageZipReader

EXAMPLES_PATH = Path(__file__).parent
RECORDING_PATH = EXAMPLES_PATH / "nanotube-example-recording.nanover.zip"
RECORDING_PATH_SWITCHING = EXAMPLES_PATH / "sim-switching-test-recording.nanover.zip"


@pytest.mark.parametrize("path,count", ((RECORDING_PATH, 1615),))
def test_n_messages(path, count):
    """
    Test examples recording have the expected number of messages.
    """
    with MessageZipReader.from_path(path) as reader:
        assert len(reader) == count
        assert sum(1 for _ in reader) == count


@pytest.mark.parametrize("path,count", ((RECORDING_PATH, 1615),))
def test_n_entries(path, count):
    entries = list(iter_recording_files(path))
    assert len(entries) == count


@pytest.mark.parametrize("path", (RECORDING_PATH, RECORDING_PATH_SWITCHING))
def test_monotonic_timestamp(path):
    """
    Test the timestamps are read correctly such that they increase monotonically.
    """
    with MessageZipReader.from_path(path) as reader:
        prev_time = 0
        for entry in reader:
            next_time = entry.metadata["timestamp"]
            assert next_time >= prev_time
            prev_time = next_time


@pytest.mark.parametrize("path", (RECORDING_PATH, RECORDING_PATH_SWITCHING))
def test_full_view_independent(path):
    """
    Test that frames and state dictionaries returned are independent and mutations don't affect subsequent frames and
    state.
    """
    KEY = "__test"

    for timestamp, frame, state in iter_full_view(path):
        assert KEY not in frame
        assert KEY not in state
        frame.values[KEY] = timestamp
        state[KEY] = timestamp


@pytest.mark.parametrize("path", (RECORDING_PATH, RECORDING_PATH_SWITCHING))
def test_full_view_no_gaps(path):
    """
    Test that frames and state dictionaries are never None and continue to persist fields present only at the start
    of the recording.
    """
    # skip past initial setup
    iterator = iter_full_view(path)
    next(iterator)
    next(iterator)

    for timestamp, frame, state in iterator:
        assert frame is not None
        assert state is not None

        # frame fields that appear during resets only
        assert frame.particle_count > 0
        assert len(frame.bond_pairs) > 0


def test_split_by_simulation_counter():
    """
    Test that splitting a known recording results in the correct number of parts with the expected length and expected
    properties e.g particle counts etc.
    """
    path = RECORDING_PATH_SWITCHING

    expected_count = 4
    expected_entry_counts = [104, 108, 127, 115]
    expected_particle_counts = [65, 173, 65, 173]

    views = split_by_simulation_counter(path)

    assert len(views) == expected_count

    for i in range(expected_count):
        view = views[i]
        timestamp, frame, state = view[-1]

        assert frame.simulation_counter == i
        assert len(view) == expected_entry_counts[i]
        assert frame.particle_count == expected_particle_counts[i]
