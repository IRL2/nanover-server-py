import pytest
import time
import threading

from osc_client import OscClient
from narupa.trajectory import FrameServer, FrameData

from pythonosc import dispatcher
from pythonosc.osc_server import ThreadingOSCUDPServer

SEND_INTERVAL = 1 / 100


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
    frame_server = FrameServer(address='localhost', port=0)
    yield frame_server
    frame_server.close()


@pytest.fixture
def osc_server():
    server = ThreadingOSCUDPServer(('localhost', 0), dispatcher.Dispatcher())
    threading.Thread(target=server.serve_forever, daemon=True).start()
    yield server
    server.shutdown()


@pytest.fixture
def server_dispatcher_client(frame_server, osc_server):
    osc_port = osc_server.socket.getsockname()[1]
    client = OscClient(osc_address='localhost', osc_port=osc_port,
                       traj_address='localhost', traj_port=frame_server.port,
                       message_generator=simple_frame_to_message,
                       send_interval=SEND_INTERVAL)
    threading.Thread(target=client.run, daemon=True).start()
    yield frame_server, osc_server, client
    client.close()


def test_test(server_dispatcher_client, simple_frame_data):
    frame_server, osc_server, osc_client = server_dispatcher_client

    test_address = "/test"
    send_message = 69
    recv_message = None

    simple_frame_data.values[test_address] = send_message

    def recv_test(address, message):
        nonlocal recv_message
        recv_message = message

    osc_server.dispatcher.map(test_address, recv_test)
    frame_server.send_frame(frame_data=simple_frame_data, frame_index=0)

    time.sleep(SEND_INTERVAL * 2)

    assert recv_message == send_message
