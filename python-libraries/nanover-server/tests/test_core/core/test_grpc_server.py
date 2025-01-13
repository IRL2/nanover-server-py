import pytest
from nanover.core import NanoverServer
from nanover.trajectory import FrameServer
from nanover.imd.imd_server import ImdServer

TEST_SERVERS = (FrameServer, ImdServer, NanoverServer)
TEST_PORTS = (54321, 60000, 62123)


@pytest.mark.parametrize("server_factory", TEST_SERVERS)
def test_address(server_factory):
    with server_factory(address="localhost", port=0) as frame_server:
        assert frame_server.address == "localhost"


@pytest.mark.parametrize("server_factory", TEST_SERVERS)
def test_any_port(server_factory):
    with server_factory(address="localhost", port=0) as frame_server:
        assert frame_server.port != 0


@pytest.mark.parametrize("server_factory", TEST_SERVERS)
def test_port_retained_after_close(server_factory):
    with server_factory(address="localhost", port=0) as frame_server:
        prev_port = frame_server.port
    assert prev_port == frame_server.port


@pytest.mark.serial
@pytest.mark.parametrize("server_factory", TEST_SERVERS)
@pytest.mark.parametrize("port", TEST_PORTS)
def test_specific_port(server_factory, port):
    with server_factory(address="localhost", port=port) as frame_server:
        assert frame_server.port == port


@pytest.mark.parametrize("server_factory", TEST_SERVERS)
def test_address(server_factory):
    with server_factory(address="localhost", port=0) as frame_server:
        assert frame_server.address == "localhost"


@pytest.mark.serial
@pytest.mark.parametrize("server_factory", TEST_SERVERS)
@pytest.mark.parametrize("port", TEST_PORTS)
def test_specific_port_in_use(server_factory, port):
    with server_factory(address="localhost", port=port) as frame_server1:
        assert frame_server1.port == port
        with pytest.raises(IOError):
            with server_factory(address="localhost", port=port):
                pass
