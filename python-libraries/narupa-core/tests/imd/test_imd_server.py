import pytest
from narupa.imd.imd_server import ImdServer

@pytest.fixture
def imd_server():
    server = ImdServer(address='localhost', port=54322)
    yield server
    server.close()

@pytest.fixture
def simple_interaction():
    pass
def test_server(imd_server):
    assert imd_server is not None

def test_public_interaction(imd_server):
    pass
    # imd_server.publish_interaction()