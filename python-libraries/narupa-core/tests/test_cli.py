import os
import signal
import time
import subprocess
import pytest


def is_process_running(process):
    return process.poll() is None


@pytest.mark.timeout(5)
@pytest.mark.skipif(os.name == 'nt', reason='Emulating ctrl-c in a subprocess does not work on Windows')
def test_run_multiplayer_server_runs():
    server_process = subprocess.Popen(["narupa-multiplayer"])
    time.sleep(1)
    assert is_process_running(server_process)
    server_process.send_signal(signal.SIGINT)
    server_process.wait()
    assert not is_process_running(server_process)
