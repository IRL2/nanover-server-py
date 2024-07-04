import time
from unittest.mock import mock_open, patch

from nanover.omni.cli import initialise_runner
from common import ARGON_XML_PATH, RECORDING_PATH_TRAJ, RECORDING_PATH_STATE


@patch("builtins.open", mock_open())
def test_record_opens_files():
    """
    Test that the expected files are opened for writing during recording.
    """
    with initialise_runner(["--record", "test", "--port", "0"]):
        pass

    open.assert_any_call("test.traj", "wb")
    open.assert_any_call("test.state", "wb")


def test_cycle_multiple_sims():
    """
    Test that multiple sims can be given and cycled through.
    """
    with initialise_runner(
        [
            "--port",
            "0",
            "--omm",
            str(ARGON_XML_PATH),
            "--ase-omm",
            str(ARGON_XML_PATH),
            "--playback",
            str(RECORDING_PATH_TRAJ),
            str(RECORDING_PATH_STATE),
        ]
    ) as runner:
        for _ in range(10):
            runner.next()
            time.sleep(0.1)
