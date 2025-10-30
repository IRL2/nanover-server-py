import pytest
import threading

from nanover.app import NanoverImdApplication
from nanover.testing import assert_equal_soon
from nanover.websocket import NanoverImdClient
from osc_client import OscClient
from nanover.trajectory import FrameData

from pythonosc import dispatcher
from pythonosc.osc_server import ThreadingOSCUDPServer

# See https://github.com/attwad/python-osc/issues/109
IPV4_LOCALHOST = "127.0.0.1"
OSC_SEND_INTERVAL = 1 / 100


def simple_frame_to_message(frame):
    try:
        yield "/test", frame["/test"]
    except KeyError:
        pass


@pytest.fixture
def simple_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data["indices"] = [0, 1, 3]
    basic_frame_data["string"] = "str"
    basic_frame_data["bool"] = False
    return basic_frame_data


@pytest.fixture
def app_server():
    """
    Provide a server hosting on an available port on localhost.
    """
    with NanoverImdApplication.basic_server(port=0) as app_server:
        yield app_server


@pytest.fixture
def osc_server():
    """
    Provide an OSC server hosting on an available port on localhost.
    """
    try:
        server = ThreadingOSCUDPServer((IPV4_LOCALHOST, 0), dispatcher.Dispatcher())
        threading.Thread(target=server.serve_forever, daemon=True).start()
        yield server
    finally:
        server.shutdown()


@pytest.fixture
def frame_osc_converter(app_server, osc_server):
    """
    Provide a frame server, OSC server, and a client that is connected to both
    of them.
    """
    osc_port = osc_server.socket.getsockname()[1]
    nanover_client = NanoverImdClient.from_app_server(app_server)
    with OscClient(
        nanover_client,
        osc_address=(IPV4_LOCALHOST, osc_port),
        message_generator=simple_frame_to_message,
        osc_send_interval=OSC_SEND_INTERVAL,
    ) as client:
        threading.Thread(target=client.run, daemon=True).start()
        yield app_server, osc_server, client


def test_transmission(frame_osc_converter, simple_frame_data):
    """
    Test that OscClient receiving frames can trigger the sending OSC messages.
    """
    app_server, osc_server, osc_client = frame_osc_converter

    test_address = "/test"
    send_message = "hello"
    recv_message = None

    simple_frame_data[test_address] = send_message

    def recv_test(address, message):
        nonlocal recv_message
        recv_message = message

    osc_server.dispatcher.map(test_address, recv_test)
    app_server.frame_publisher.send_frame(frame=simple_frame_data)

    assert_equal_soon(
        lambda: recv_message,
        lambda: send_message,
    )
