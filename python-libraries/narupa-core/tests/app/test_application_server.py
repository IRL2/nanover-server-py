import pytest
from narupa.app import NarupaApplicationServer


@pytest.mark.serial
def test_run_two_servers_default_port():
    with NarupaApplicationServer.basic_server():
        with pytest.raises(IOError):
            with NarupaApplicationServer.basic_server():
                pass


def test_run_two_servers_same_port():
    with NarupaApplicationServer.basic_server(port=0) as server:
        with pytest.raises(IOError):
            with NarupaApplicationServer.basic_server(port=server.port):
                pass
