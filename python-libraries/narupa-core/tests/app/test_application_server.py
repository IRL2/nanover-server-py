import pytest
from nanover.app import NanoVerApplicationServer


@pytest.mark.serial
def test_run_two_servers_default_port():
    with NanoVerApplicationServer.basic_server():
        with pytest.raises(IOError):
            with NanoVerApplicationServer.basic_server():
                pass


def test_run_two_servers_same_port():
    with NanoVerApplicationServer.basic_server(port=0) as server:
        with pytest.raises(IOError):
            with NanoVerApplicationServer.basic_server(port=server.port):
                pass
