import os
import time

import pytest

from nanover.app.cli.omni_cli import initialise_runner, handle_user_arguments
from common import ARGON_XML_PATH, RECORDING_PATH


@pytest.mark.serial
def test_record_writes_files(tmp_path):
    """
    Test that the expected files created during recording.
    """
    arguments = handle_user_arguments(["--record", str(tmp_path / "test")])

    with initialise_runner(arguments):
        pass

    assert set(os.listdir(tmp_path)) == {"test.nanover.zip"}


@pytest.mark.serial
def test_cycle_multiple_sims():
    """
    Test that multiple sims can be given and cycled through.
    """
    arguments = handle_user_arguments(
        [
            "--omm",
            str(ARGON_XML_PATH),
            "--playback",
            str(RECORDING_PATH),
        ]
    )

    with initialise_runner(arguments) as runner:
        for _ in range(4):
            runner.next()
            time.sleep(0.05)
