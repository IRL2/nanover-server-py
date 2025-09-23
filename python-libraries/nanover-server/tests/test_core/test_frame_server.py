import pytest

from nanover.trajectory.frame_data2 import FrameData


@pytest.fixture
def simple_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data["indices"] = [0, 1, 3]
    basic_frame_data["string"] = "str"
    basic_frame_data["bool"] = False
    return basic_frame_data


@pytest.fixture
def disjoint_frame_data():
    basic_frame_data = FrameData()
    basic_frame_data["strings"] = ["a", "b", "d"]
    basic_frame_data["number"] = 16.5
    return basic_frame_data
