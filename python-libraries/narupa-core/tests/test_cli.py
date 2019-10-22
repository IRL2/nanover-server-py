import signal
import time
import subprocess
import pytest


def is_process_running(process):
    return process.poll() is None


@pytest.mark.timeout(5)
def test_run_multiplayer_server_runs():
    server_process = subprocess.Popen(["narupa-multiplayer"])
    time.sleep(1)
    assert is_process_running(server_process)
    server_process.send_signal(signal.CTRL_C_EVENT)
    server_process.wait()



