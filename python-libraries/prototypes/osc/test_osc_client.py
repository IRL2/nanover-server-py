import pytest
import time
import threading

from osc_client import OscClient
from narupa.trajectory import FrameServer, FrameData
from narupa.app.client import DEFAULT_SUBSCRIPTION_INTERVAL

from pythonosc import dispatcher
from pythonosc.osc_server import ThreadingOSCUDPServer

OSC_SEND_INTERVAL = 1 / 100


def simple_frame_to_message(frame):
    yield "/test", frame.values["/test"]


@pytest.fixture
def simple_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data.arrays["indices"] = [0, 1, 3]
    basic_frame_data.values["string"] = "str"
    basic_frame_data.values["bool"] = False
    return basic_frame_data


@pytest.fixture
def frame_server():
    """
    Provide a frame server hosting on an available port on localhost.
    """
    with FrameServer(address='localhost', port=0) as frame_server:
        yield frame_server


@pytest.fixture
def osc_server():
    """
    Provide an OSC server hosting on an available port on localhost.
    """
    try:
        server = ThreadingOSCUDPServer(('localhost', 0), dispatcher.Dispatcher())
        threading.Thread(target=server.serve_forever, daemon=True).start()
        yield server
    finally:
        server.shutdown()


@pytest.fixture
def frame_osc_converter(frame_server, osc_server):
    """
    Provide a frame server, OSC server, and a client that is connected to both
    of them.
    """
    osc_port = osc_server.socket.getsockname()[1]
    with OscClient(osc_address='localhost', osc_port=osc_port,
                   traj_address='localhost', traj_port=frame_server.port,
                   message_generator=simple_frame_to_message,
                   send_interval=OSC_SEND_INTERVAL) as client:
        threading.Thread(target=client.run, daemon=True).start()
        yield frame_server, osc_server, client


def test_transmission(frame_osc_converter, simple_frame_data):
    """
    Test that OscClient receiving frames can trigger the sending OSC messages.
    """
    frame_server, osc_server, osc_client = frame_osc_converter

    test_address = "/test"
    send_message = "hello"
    recv_message = None

    simple_frame_data.values[test_address] = send_message

    def recv_test(address, message):
        nonlocal recv_message
        recv_message = message

    osc_server.dispatcher.map(test_address, recv_test)
    frame_server.send_frame(frame_data=simple_frame_data, frame_index=0)

    time.sleep(DEFAULT_SUBSCRIPTION_INTERVAL * 2)
    time.sleep(OSC_SEND_INTERVAL * 2)

    assert recv_message == send_message
