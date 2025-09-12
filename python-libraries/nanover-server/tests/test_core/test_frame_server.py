from typing import Iterable
from contextlib import contextmanager

import numpy
import pytest
import time

from grpc import RpcError, StatusCode

from nanover.app import NanoverImdApplication
from nanover.trajectory import FrameData
from nanover.trajectory.frame_data import (
    SERVER_TIMESTAMP,
    SIMULATION_COUNTER,
    FRAME_INDEX,
)
from numpy import average

from nanover.websocket.client.app_client import NanoverImdClient
from .utilities.test_timing import TIMING_TOLERANCE, COMMON_INTERVALS

FRAME_DATA_VARIABLE_KEYS = (SERVER_TIMESTAMP, SIMULATION_COUNTER, FRAME_INDEX)
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
