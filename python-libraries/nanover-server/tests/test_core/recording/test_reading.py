from pathlib import Path

import pytest

from nanover.recording.reading import (
    iter_trajectory_file,
    iter_state_file,
    iter_recording_files,
    iter_full_view,
    split_by_simulation_counter,
    MessageRecordingReader,
)

EXAMPLES_PATH = Path(__file__).parent
RECORDING_PATH_TRAJ = EXAMPLES_PATH / "nanotube-example-recording.traj"
RECORDING_PATH_STATE = EXAMPLES_PATH / "nanotube-example-recording.state"

RECORDING_PATH_TRAJ_SWITCHING = EXAMPLES_PATH / "sim-switching-test-recording.traj"
RECORDING_PATH_STATE_SWITCHING = EXAMPLES_PATH / "sim-switching-test-recording.state"


@pytest.mark.parametrize(
    "path,count", ((RECORDING_PATH_TRAJ, 930), (RECORDING_PATH_STATE, 685))
)
def test_n_messages(path, count):
    """
    Test examples recording have the expected number of messages.
    """
    with MessageRecordingReader.from_path(path) as reader:
        assert len(reader) == count
        assert sum(1 for _ in reader) == count


@pytest.mark.parametrize(
    "traj_path,state_path,count", ((RECORDING_PATH_TRAJ, RECORDING_PATH_STATE, 1615),)
)
def test_n_entries(traj_path, state_path, count):
    entries = list(iter_recording_files(traj=traj_path, state=state_path))
    assert len(entries) == count


@pytest.mark.parametrize("path", (RECORDING_PATH_TRAJ, RECORDING_PATH_STATE))
def test_monotonic_timestamp(path):
    """
    Test the timestamps are read correctly such that they increase monotonically.
    """
    with MessageRecordingReader.from_path(path) as reader:
        prev_time = 0
        for entry in reader:
            assert entry.timestamp >= prev_time
            prev_time = entry.timestamp


@pytest.mark.parametrize(
    "traj_path,state_path", ((RECORDING_PATH_TRAJ, RECORDING_PATH_STATE),)
)
def test_full_view_independent(traj_path, state_path):
    """
    Test that frames and state dictionaries returned are independent and mutations don't affect subsequent frames and
    state.
    """
    KEY = "__test"

    for timestamp, frame, state in iter_full_view(traj=traj_path, state=state_path):
        assert KEY not in frame
        assert KEY not in state
        frame.values[KEY] = timestamp
        state[KEY] = timestamp


@pytest.mark.parametrize(
    "traj_path,state_path", ((RECORDING_PATH_TRAJ, RECORDING_PATH_STATE),)
)
def test_full_view_no_gaps(traj_path, state_path):
    """
    Test that frames and state dictionaries are never None and continue to persist fields present only at the start
    of the recording.
    """
    # skip past initial setup
    iterator = iter_full_view(traj=traj_path, state=state_path)
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
    traj_path = RECORDING_PATH_TRAJ_SWITCHING
    state_path = RECORDING_PATH_STATE_SWITCHING

    expected_count = 4
    expected_entry_counts = [104, 108, 127, 115]
    expected_particle_counts = [65, 173, 65, 173]

    views = split_by_simulation_counter(traj=traj_path, state=state_path)

    assert len(views) == expected_count

    for i in range(expected_count):
        view = views[i]
        timestamp, frame, state = view[-1]

        assert frame.simulation_counter == i
        assert len(view) == expected_entry_counts[i]
        assert frame.particle_count == expected_particle_counts[i]
