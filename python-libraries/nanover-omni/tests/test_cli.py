import os
import time

import pytest

from nanover.omni.cli import initialise_runner, handle_user_arguments
from common import ARGON_XML_PATH, RECORDING_PATH_TRAJ, RECORDING_PATH_STATE


@pytest.mark.serial
def test_record_writes_files(tmp_path):
    """
    Test that the expected files created during recording.
    """
    arguments = handle_user_arguments(
        ["--record", str(tmp_path / "test"), "--port", "0"]
    )

    with initialise_runner(arguments):
        pass

    assert set(os.listdir(tmp_path)) == {"test.traj", "test.state"}


@pytest.mark.serial
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
            "--playback",
            str(RECORDING_PATH_TRAJ),
            str(RECORDING_PATH_STATE),
        ]
    )

    with initialise_runner(arguments) as runner:
        for _ in range(4):
            runner.next()
            time.sleep(0.05)
