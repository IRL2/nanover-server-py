import os
import time
import subprocess
import pytest


def is_process_running(process):
    return process.poll() is None


@pytest.mark.timeout(5)
def test_run_multiplayer_server_runs():
    # assumes test is colocated with the example which is not ideal but there's
    # no easy alternative for now
    root = os.path.dirname(os.path.realpath(__file__))
    path = os.path.join(root, "run_multiplayer_server.py")
    server_process = subprocess.Popen(["python", path])

    time.sleep(1)
    assert is_process_running(server_process)
    server_process.terminate()
