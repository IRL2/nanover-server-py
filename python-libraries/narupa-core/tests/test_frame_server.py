from unittest.mock import Mock

import pytest
import time

from narupa.trajectory import FrameServer, FrameClient, FrameData


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
    server = FrameServer(address='localhost', port=54321)
    yield server
    server.close()


@pytest.fixture
def frame_client():
    client = FrameClient(address='localhost', port=54321)
    yield client
    client.close()


def test_blankdata_lateclient(frame_server, frame_client):
    mock = Mock()

    frame_server.send_frame(0, FrameData())

    frame_client.subscribe_frames_async(mock.callback)

    time.sleep(0.5)

    mock.callback.assert_called_once()


def test_blankdata_earlyclient(frame_server, frame_client):
    mock = Mock()

    frame_client.subscribe_frames_async(mock.callback)

    frame_server.send_frame(0, FrameData())

    time.sleep(0.5)

    mock.callback.assert_called_once()


# Checks the transmitted data is correct
def test_data_earlyclient(frame_server, frame_client, simple_frame_data):
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    frame_client.subscribe_frames_async(callback)
    # It takes time to actually subscribe. During that time, the server can
    # already have yielded the frame from the next instruction. We therefore
    # need to wait for the subscription go go through before we send a frame.
    time.sleep(0.1)

    frame_server.send_frame(0, simple_frame_data)

    time.sleep(0.5)

    assert result == simple_frame_data


def test_data_lateclient(frame_server, frame_client, simple_frame_data):
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    frame_server.send_frame(0, simple_frame_data)

    frame_client.subscribe_frames_async(callback)

    time.sleep(0.5)

    assert result == simple_frame_data


def test_data_disjoint(frame_server, frame_client, simple_frame_data, disjoint_frame_data,
                       simple_and_disjoint_frame_data):
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    frame_server.send_frame(0, simple_frame_data)
    frame_server.send_frame(1, disjoint_frame_data)

    frame_client.subscribe_frames_async(callback)

    time.sleep(0.5)

    assert result == simple_and_disjoint_frame_data


def test_data_overlap(frame_server, frame_client, simple_frame_data, overlap_frame_data, simple_and_overlap_frame_data):
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    frame_server.send_frame(0, simple_frame_data)
    frame_server.send_frame(1, overlap_frame_data)

    frame_client.subscribe_frames_async(callback)

    time.sleep(0.5)

    assert result == simple_and_overlap_frame_data


def test_slow_frame_publishing(frame_server, frame_client, simple_frame_data):
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    future = frame_client.subscribe_frames_async(callback)
    time.sleep(0.1)

    for i in range(5):
        time.sleep(0.5)
        frame_server.send_frame(i, simple_frame_data)

    time.sleep(0.5)
    future.cancel()
    # this result would throw an exception if an exception was raised during subscription.
    subscribe_result = future.result()
    assert subscribe_result is None
    assert result == simple_frame_data
