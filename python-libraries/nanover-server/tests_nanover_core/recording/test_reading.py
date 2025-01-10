from pathlib import Path

import pytest

from nanover.recording.reading import (
    iter_trajectory_file,
    iter_state_file,
    iter_recording_buffers,
    iter_recording_files,
)

EXAMPLES_PATH = Path(__file__).parent
RECORDING_PATH_TRAJ = EXAMPLES_PATH / "nanotube-example-recording.traj"
RECORDING_PATH_STATE = EXAMPLES_PATH / "nanotube-example-recording.state"


@pytest.mark.parametrize("path,count", ((RECORDING_PATH_TRAJ, 930),))
def test_n_frames(path, count):
    """
    Test an example recording has the expected number of frames.
    """
    frames = list(iter_trajectory_file(path))
    assert len(frames) == count


@pytest.mark.parametrize("path,count", ((RECORDING_PATH_STATE, 685),))
def test_n_updates(path, count):
    """
    Test an example recording has the expected number of updates.
    """
    updates = list(iter_state_file(path))
    assert len(updates) == count


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
    with open(path, "rb") as infile:
        prev_time = 0
        for next_time, _ in iter_recording_buffers(infile):
            assert next_time >= prev_time
            prev_time = next_time
