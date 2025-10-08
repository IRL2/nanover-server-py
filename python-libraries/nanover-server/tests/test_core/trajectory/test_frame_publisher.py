"""
Unit tests for `nanover.trajectory.frame_publisher
"""

from concurrent.futures import ThreadPoolExecutor
import itertools
import pytest

from nanover.trajectory import FramePublisher, FrameData
from nanover.trajectory.keys import SERVER_TIMESTAMP
from nanover.utilities.cli import CancellationToken


def test_send_frame_data():
    publisher = FramePublisher()
    publisher.send_frame(3, FrameData())
    assert publisher.last_frame.frame_index == 3
    assert SERVER_TIMESTAMP in publisher.last_frame


@pytest.mark.timeout(1)
def test_cancellation_ends_empty_stream():
    """
    Test that a cancelled a frame stream yields no frames and terminates.
    """
    publisher = FramePublisher()
    cancellation = CancellationToken()

    stream = publisher.subscribe_latest_frames(
        frame_interval=0,
        cancellation=cancellation,
    )

    cancellation.cancel()

    assert sum(1 for _ in stream) == 0


@pytest.mark.timeout(1)
@pytest.mark.parametrize("count", list(range(5)))
def test_cancellation_ends_stream_immediately(count):
    """
    Test that a cancelled a frame stream yields no frames and terminates, even though frames were sent.
    """
    publisher = FramePublisher()
    cancellation = CancellationToken()

    stream = publisher.subscribe_latest_frames(
        frame_interval=0,
        cancellation=cancellation,
    )

    for i in range(count):
        publisher.send_frame(i, FrameData())

    cancellation.cancel()

    assert sum(1 for _ in stream) == 0


@pytest.mark.timeout(1)
def test_cancellation_ends_started_stream():
    """
    Test that a cancelled a frame stream yields no more frames and terminates, even after consuming some frames and then
    sending more.
    """
    count = 5

    publisher = FramePublisher()
    cancellation = CancellationToken()

    stream = publisher.subscribe_latest_frames(
        frame_interval=0,
        cancellation=cancellation,
    )

    for i in range(count):
        publisher.send_frame(i, FrameData())
        next(stream)

    for _ in range(count):
        publisher.send_frame(i, FrameData())

    cancellation.cancel()

    assert sum(1 for _ in stream) == 0
