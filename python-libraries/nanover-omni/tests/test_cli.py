import time
from unittest.mock import mock_open, patch

from nanover.omni.cli import initialise_runner, handle_user_arguments
from common import ARGON_XML_PATH, RECORDING_PATH_TRAJ, RECORDING_PATH_STATE


@patch("builtins.open", mock_open())
def test_record_opens_files():
    """
    Test that the expected files are opened for writing during recording.
    """
    arguments = handle_user_arguments(["--record", "test", "--port", "0"])

    with initialise_runner(arguments):
        pass

    open.assert_any_call("test.traj", "wb")
    open.assert_any_call("test.state", "wb")


def test_cycle_multiple_sims():
    """
    Test that multiple sims can be given and cycled through.
    """
    arguments = handle_user_arguments(
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
    )

    with initialise_runner(arguments) as runner:
        for _ in range(4):
            runner.next()
            time.sleep(0.05)
