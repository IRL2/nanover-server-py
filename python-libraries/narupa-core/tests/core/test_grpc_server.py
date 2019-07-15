import pytest
from narupa.trajectory import FrameServer
from narupa.imd.imd_server import ImdServer
from narupa.multiplayer.multiplayer_server import MultiplayerServer

TEST_SERVERS = (FrameServer, ImdServer, MultiplayerServer)
TEST_PORTS = (54321, 60000, 62123)


@pytest.mark.parametrize('server_factory', TEST_SERVERS)
def test_any_port(server_factory):
    frame_server = server_factory(address='localhost', port=0)
    try:
        assert frame_server.port != 0
    finally:
        frame_server.close()


@pytest.mark.parametrize('server_factory', TEST_SERVERS)
def test_port_retained_after_close(server_factory):
    frame_server = server_factory(address='localhost', port=0)
    prev_port = frame_server.port
    frame_server.close()
    assert prev_port == frame_server.port


@pytest.mark.parametrize('server_factory', TEST_SERVERS)
@pytest.mark.parametrize('port', TEST_PORTS)
def test_specific_port(server_factory, port):
    frame_server = server_factory(address='localhost', port=port)
    try:
        assert frame_server.port == port
    finally:
        frame_server.close()


@pytest.mark.parametrize('server_factory', TEST_SERVERS)
@pytest.mark.parametrize('port', TEST_PORTS)
def test_specific_port_in_use(server_factory, port):
    frame_server1 = server_factory(address='localhost', port=port)
    frame_server2 = server_factory(address='localhost', port=port)
    try:
        assert frame_server1.port == port and frame_server2.port == 0
    finally:
        frame_server1.close()
        frame_server2.close()



