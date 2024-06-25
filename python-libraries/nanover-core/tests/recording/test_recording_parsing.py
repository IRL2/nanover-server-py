from pathlib import Path

from nanover.recording.parsing import iter_trajectory_file, iter_state_file

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
