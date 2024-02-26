from typing import Iterable
from contextlib import contextmanager
from unittest.mock import Mock

import pytest
import time

from grpc import RpcError, StatusCode
from nanover.trajectory import FrameServer, FrameClient, FrameData
from nanover.trajectory.frame_data import SERVER_TIMESTAMP
from numpy import average

SUBSCRIBE_METHODS = ("subscribe_frames_async", "subscribe_last_frames_async")
FRAME_DATA_VARIABLE_KEYS = (SERVER_TIMESTAMP,)
IMMEDIATE_REPLY_WAIT_TIME = 0.01


def assert_framedata_equal(
    left: FrameData,
    right: FrameData,
    ignore_keys: Iterable[str] = FRAME_DATA_VARIABLE_KEYS,
):
    """
    Raise an :exc:`AssertError` if the two frames are not equal.

    One can ignore keys from the comparison by listing them in the `ignore_key`
    argument.
    """
    left = left.copy()
    right = right.copy()

    for key in ignore_keys:
        del left[key], right[key]

    assert left == right


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
    basic_frame_data.arrays["strings"] = ["a", "b", "d"]
    basic_frame_data.values["number"] = 16.5
    return basic_frame_data


@pytest.fixture
def overlap_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data.arrays["strings"] = ["a", "b", "d"]
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
    basic_frame_data.arrays["strings"] = ["a", "b", "d"]
    basic_frame_data.values["number"] = 16.5
    return basic_frame_data


@pytest.fixture
def simple_and_overlap_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data.values["string"] = "str"
    basic_frame_data.arrays["strings"] = ["a", "b", "d"]
    basic_frame_data.arrays["indices"] = [6, 8, 11]
    basic_frame_data.values["number"] = 16.5
    basic_frame_data.values["bool"] = True
    return basic_frame_data


@pytest.fixture
def frame_server():
    with FrameServer(address="localhost", port=0) as frame_server:
        yield frame_server


@pytest.fixture
def frame_server_client_pair(frame_server):
    client = FrameClient.insecure_channel(address="localhost", port=frame_server.port)
    yield frame_server, client
    client.close()


@pytest.mark.parametrize("subscribe_method", SUBSCRIBE_METHODS)
def test_blankdata_lateclient(frame_server_client_pair, subscribe_method):
    frame_server, frame_client = frame_server_client_pair

    mock = Mock()

    frame_server.send_frame(0, FrameData())

    getattr(frame_client, subscribe_method)(mock.callback)

    time.sleep(0.1)

    mock.callback.assert_called_once()


@pytest.mark.parametrize("subscribe_method", SUBSCRIBE_METHODS)
def test_blankdata_earlyclient(frame_server_client_pair, subscribe_method):
    frame_server, frame_client = frame_server_client_pair

    mock = Mock()

    getattr(frame_client, subscribe_method)(mock.callback)

    frame_server.send_frame(0, FrameData())

    time.sleep(0.1)

    mock.callback.assert_called_once()


# Checks the transmitted data is correct
@pytest.mark.parametrize("subscribe_method", SUBSCRIBE_METHODS)
def test_data_earlyclient(
    frame_server_client_pair, simple_frame_data, subscribe_method
):
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

    time.sleep(0.1)

    assert_framedata_equal(result, simple_frame_data)


@pytest.mark.parametrize("subscribe_method", SUBSCRIBE_METHODS)
def test_data_lateclient(frame_server_client_pair, simple_frame_data, subscribe_method):
    frame_server, frame_client = frame_server_client_pair
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    frame_server.send_frame(0, simple_frame_data)

    getattr(frame_client, subscribe_method)(callback)

    time.sleep(0.1)
    assert SERVER_TIMESTAMP in result.values

    assert_framedata_equal(result, simple_frame_data)


@pytest.mark.parametrize("subscribe_method", SUBSCRIBE_METHODS)
def test_data_disjoint(
    frame_server_client_pair,
    simple_frame_data,
    disjoint_frame_data,
    simple_and_disjoint_frame_data,
    subscribe_method,
):
    frame_server, frame_client = frame_server_client_pair
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    frame_server.send_frame(0, simple_frame_data)
    frame_server.send_frame(1, disjoint_frame_data)

    getattr(frame_client, subscribe_method)(callback)

    time.sleep(0.1)
    assert SERVER_TIMESTAMP in result.values

    assert_framedata_equal(result, simple_and_disjoint_frame_data)


@pytest.mark.parametrize("subscribe_method", SUBSCRIBE_METHODS)
def test_data_overlap(
    frame_server_client_pair,
    simple_frame_data,
    overlap_frame_data,
    simple_and_overlap_frame_data,
    subscribe_method,
):
    frame_server, frame_client = frame_server_client_pair
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    frame_server.send_frame(0, simple_frame_data)
    frame_server.send_frame(1, overlap_frame_data)

    getattr(frame_client, subscribe_method)(callback)

    time.sleep(0.1)
    assert SERVER_TIMESTAMP in result.values

    assert_framedata_equal(result, simple_and_overlap_frame_data)


@contextmanager
def raises_rpc_cancelled():
    """
    Silently ignore an RpcError exception with the CANCELLED status code.
    """
    try:
        yield
    except RpcError as e:
        if e._state.code != StatusCode.CANCELLED:
            raise e


@pytest.mark.parametrize("subscribe_method", SUBSCRIBE_METHODS)
def test_slow_frame_publishing(
    frame_server_client_pair, simple_frame_data, subscribe_method
):
    frame_server, frame_client = frame_server_client_pair
    result = None

    def callback(frame, **kwargs):
        nonlocal result
        result = frame

    future = getattr(frame_client, subscribe_method)(callback)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    for i in range(5):
        time.sleep(0.1)
        frame_server.send_frame(i, simple_frame_data)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    # TODO there is no way to cancel the subscription stream...
    frame_client.close()

    with raises_rpc_cancelled():
        future.result()

    assert_framedata_equal(result, simple_frame_data)


def test_subscribe_latest_frames_sends_latest_frame(
    frame_server_client_pair, simple_frame_data
):
    frame_server, frame_client = frame_server_client_pair

    frame_interval = 1 / 30
    first_index = None

    def callback(frame, frame_index):
        nonlocal first_index
        if first_index is None:
            first_index = frame_index

    frame_client.subscribe_last_frames_async(callback)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    for i in range(5):
        frame_server.send_frame(i, simple_frame_data)

    time.sleep(2 * frame_interval)
    assert first_index == 4


def test_subscribe_latest_frames_aggregates_frames(
    frame_server_client_pair,
    simple_frame_data,
):
    """
    Test that data that exists only in intermediate frames (those frames that
    are never "latest" at the time of sending) is aggregated with the latest
    frame instead of being lost.
    """
    frame_server, frame_client = frame_server_client_pair

    frame_interval = 1 / 30
    latest_frame = None
    latest_index = None

    initial = FrameData()
    initial.values["every.frame"] = "initial"
    initial.values["some.frames"] = "initial"

    intermediate = FrameData()
    intermediate.values["every.frame"] = "intermediate"
    intermediate.values["some.frames"] = "intermediate"

    latest = FrameData()
    latest.values["every.frame"] = "latest"

    def callback(frame, frame_index):
        nonlocal latest_frame, latest_index
        latest_frame = frame
        latest_index = frame_index

    frame_server.send_frame(0, initial)

    frame_client.subscribe_last_frames_async(callback)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    # add two frames at once, because of the subscription interval only one
    # frame will be sent, and we expect it to be the aggregate of both frames
    frame_server.send_frame(1, intermediate)
    frame_server.send_frame(2, latest)

    time.sleep(2 * frame_interval)
    assert latest_frame.values["some.frames"] == intermediate.values["some.frames"]
    assert latest_frame.values["every.frame"] == latest.values["every.frame"]
    assert latest_index == 2


@pytest.mark.parametrize("subscribe_method", SUBSCRIBE_METHODS)
@pytest.mark.parametrize("frame_interval", (1 / 10, 1 / 30, 1 / 60))
def test_subscribe_frames_frame_interval(
    frame_server_client_pair, simple_frame_data, subscribe_method, frame_interval
):
    """
    Test that when using frame subscription methods with an interval, frames are
    sent, on average, at the specified rate.
    """
    frame_server, frame_client = frame_server_client_pair
    subscribe_frames_async = getattr(frame_client, subscribe_method)

    tolerance = 5 / 1000  # 5 milliseconds
    frame_limit = 30
    time_limit = frame_limit * frame_interval * 1.5
    receive_times = []

    def callback(frame, frame_index):
        receive_times.append(time.monotonic())

        if frame_index < frame_limit:
            frame_server.send_frame(frame_index + 1, simple_frame_data)

    subscribe_frames_async(callback, frame_interval)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    frame_server.send_frame(0, simple_frame_data)

    time.sleep(time_limit)
    intervals = [
        receive_times[i + 1] - receive_times[i] for i in range(frame_limit - 1)
    ]
    assert abs(average(intervals) - frame_interval) < tolerance
