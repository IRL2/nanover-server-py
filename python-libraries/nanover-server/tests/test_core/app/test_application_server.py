import pytest
from nanover.app import NanoverApplicationServer


@pytest.mark.serial
def test_run_two_servers_default_port():
    with NanoverApplicationServer.basic_server():
        with pytest.raises(IOError):
            with NanoverApplicationServer.basic_server():
                pass


def test_run_two_servers_same_port():
    with NanoverApplicationServer.basic_server(port=0) as server:
        with pytest.raises(IOError):
            with NanoverApplicationServer.basic_server(port=server.port):
                pass
