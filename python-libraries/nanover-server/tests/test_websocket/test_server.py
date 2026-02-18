from pytest import raises

from nanover.app import NanoverImdApplication
from nanover.websocket.server import WebSocketServer, _get_server_port


def test_port_in_use():
    """
    Test that attempting to reuse a port will raise an exception including the port number in question.
    """
    with WebSocketServer(NanoverImdApplication()) as server:
        first = server.serve()
        port = _get_server_port(first)
        with raises(IOError, match=f"{port}"):
            server.serve(port=port)
