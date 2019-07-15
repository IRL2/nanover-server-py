from unittest.mock import Mock

import pytest
import time

from grpc import RpcError, StatusCode
from narupa.trajectory import FrameServer, FrameClient, FrameData

SUBSCRIBE_METHODS = ('subscribe_frames_async', 'subscribe_last_frames_async')

@pytest.fixture
def simple_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data.arrays["indices"] = [0, 1, 3]
    basic_frame_data.values["string"] = "str"
    basic_frame_data.values["bool"] = False
    return basic_frame_data


@pytest.fixture
def disjoint_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data.arrays["strings"] = ['a', 'b', 'd']
    basic_frame_data.values["number"] = 16.5
    return basic_frame_data


@pytest.fixture
def overlap_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data.arrays["strings"] = ['a', 'b', 'd']
    basic_frame_data.arrays["indices"] = [6, 8, 11]
    basic_frame_data.values["number"] = 16.5
    basic_frame_data.values["bool"] = True
    return basic_frame_data


@pytest.fixture
def simple_and_disjoint_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data.arrays["indices"] = [0, 1, 3]
    basic_frame_data.values["string"] = "str"
    basic_frame_data.values["bool"] = False
    basic_frame_data.arrays["strings"] = ['a', 'b', 'd']
    basic_frame_data.values["number"] = 16.5
    return basic_frame_data


@pytest.fixture
def simple_and_overlap_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data.values["string"] = "str"
    basic_frame_data.arrays["strings"] = ['a', 'b', 'd']
    basic_frame_data.arrays["indices"] = [6, 8, 11]
    basic_frame_data.values["number"] = 16.5
    basic_frame_data.values["bool"] = True
    return basic_frame_data

@pytest.fixture
def frame_server():
    frame_server = FrameServer(address='localhost', port=0)
    yield frame_server
    frame_server.close()

@pytest.fixture
def frame_server_client_pair(frame_server):
    client = FrameClient(address='localhost', port=frame_server.port)
    yield frame_server, client
    client.close()

@pytest.mark.parametrize('subscribe_method', SUBSCRIBE_METHODS)
def test_blankdata_lateclient(frame_server_client_pair, subscribe_method):
    frame_server, frame_client = frame_server_client_pair
    
    mock = Mock()

    frame_server.send_frame(0, FrameData())

    getattr(frame_client, subscribe_method)(mock.callback)

    time.sleep(0.5)

    mock.callback.assert_called_once()


@pytest.mark.parametrize('subscribe_method', SUBSCRIBE_METHODS)
def test_blankdata_earlyclient(frame_server_client_pair, subscribe_method):
    frame_server, frame_client = frame_server_client_pair
    
    mock = Mock()

    getattr(frame_client, subscribe_method)(mock.callback)

    frame_server.send_frame(0, FrameData())

    time.sleep(0.5)

    mock.callback.assert_called_once()


# Checks the transmitted data is correct
@pytest.mark.parametrize('subscribe_method', SUBSCRIBE_METHODS)
def test_data_earlyclient(frame_server_client_pair, simple_frame_data,
                          subscribe_method):
    frame_server, frame_client = frame_server_client_pair
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    getattr(frame_client, subscribe_method)(callback)
    # It takes time to actually subscribe. During that time, the server can
    # already have yielded the frame from the next instruction. We therefore
    # need to wait for the subscription go go through before we send a frame.
    time.sleep(0.1)

    frame_server.send_frame(0, simple_frame_data)

    time.sleep(0.5)

    assert result == simple_frame_data


@pytest.mark.parametrize('subscribe_method', SUBSCRIBE_METHODS)
def test_data_lateclient(frame_server_client_pair, simple_frame_data,
                         subscribe_method):
    frame_server, frame_client = frame_server_client_pair
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    frame_server.send_frame(0, simple_frame_data)

    getattr(frame_client, subscribe_method)(callback)

    time.sleep(0.5)

    assert result == simple_frame_data


@pytest.mark.parametrize('subscribe_method', SUBSCRIBE_METHODS)
def test_data_disjoint(frame_server_client_pair, simple_frame_data,
                       disjoint_frame_data, simple_and_disjoint_frame_data,
                       subscribe_method):
    frame_server, frame_client = frame_server_client_pair
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    frame_server.send_frame(0, simple_frame_data)
    frame_server.send_frame(1, disjoint_frame_data)

    getattr(frame_client, subscribe_method)(callback)

    time.sleep(0.5)

    assert result == simple_and_disjoint_frame_data


@pytest.mark.parametrize('subscribe_method', SUBSCRIBE_METHODS)
def test_data_overlap(frame_server_client_pair, simple_frame_data,
                      overlap_frame_data, simple_and_overlap_frame_data,
                      subscribe_method):
    frame_server, frame_client = frame_server_client_pair
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    frame_server.send_frame(0, simple_frame_data)
    frame_server.send_frame(1, overlap_frame_data)

    getattr(frame_client, subscribe_method)(callback)

    time.sleep(0.5)

    assert result == simple_and_overlap_frame_data


@pytest.mark.parametrize('subscribe_method', SUBSCRIBE_METHODS)
def test_slow_frame_publishing(frame_server_client_pair, simple_frame_data,
                               subscribe_method):
    frame_server, frame_client = frame_server_client_pair
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    future = getattr(frame_client, subscribe_method)(callback)
    time.sleep(0.1)

    for i in range(5):
        time.sleep(0.5)
        frame_server.send_frame(i, simple_frame_data)

    time.sleep(0.5)
    # TODO there is no way to cancel the subscription stream...
    frame_client.close()
    # this result would throw an exception if an exception was raised during subscription.
    # here, we handle the expected cancelled exception, anything else is a bug.
    try:
        subscribe_result = future.result()
    except RpcError as e:
        if not e._state.code == StatusCode.CANCELLED:
            raise e

    assert result == simple_frame_data
