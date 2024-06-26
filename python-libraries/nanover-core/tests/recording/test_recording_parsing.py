from pathlib import Path

import pytest

from nanover.recording.parsing import (
    iter_trajectory_file,
    iter_state_file,
    iter_recording_buffers,
)
from nanover.recording.unpacker import Unpacker

EXAMPLES_PATH = Path(__file__).parent
RECORDING_PATH_TRAJ = EXAMPLES_PATH / "nanotube-example-recording.traj"
RECORDING_PATH_STATE = EXAMPLES_PATH / "nanotube-example-recording.state"


def test_parse_traj():
    """
    Test an example recording has the expected number of frames.
    """
    frames = list(iter_trajectory_file(RECORDING_PATH_TRAJ))
    assert len(frames) == 929


def test_parse_state():
    """
    Test an example recording has the expected number of updates.
    """
    updates = list(iter_state_file(RECORDING_PATH_STATE))
    assert len(updates) == 684


@pytest.mark.parametrize("path", (RECORDING_PATH_TRAJ, RECORDING_PATH_STATE))
def test_monotonic_timestamp(path):
    """
    Test the timestamps are read correctly such that they increase monotonically.
    """
    with open(path, "rb") as infile:
        data = infile.read()
    unpacker = Unpacker(data)

    prev_time = 0
    for next_time, _ in iter_recording_buffers(unpacker):
        assert next_time >= prev_time
        prev_time = next_time
