"""
Unit tests for `nanover.trajectory.frame_publisher
"""

from concurrent.futures import ThreadPoolExecutor
import itertools
import pytest

from nanover.trajectory import FramePublisher, FrameData
from nanover.trajectory.keys import SERVER_TIMESTAMP
from nanover.utilities.cli import CancellationToken


def test_user_queue():
    """
    The `_user_queue` works as expected for a context manager.
    """
    publisher = FramePublisher()
    assert not publisher.frame_queues.queues
    with publisher.frame_queues.one_queue(0):
        assert list(publisher.frame_queues.queues.keys()) == [0]
    assert not publisher.frame_queues.queues


def test_send_frame_data():
    publisher = FramePublisher()
    publisher.send_frame(3, FrameData())
    assert publisher.last_frame.frame_index == 3
    assert SERVER_TIMESTAMP in publisher.last_frame


def test_get_new_request_id_serial():
    """
    `_get_new_request_id` works in serial.
    """
    number_of_ids = 5
    publisher = FramePublisher()
    obtained = [publisher._get_new_request_id() for _ in range(number_of_ids)]
    assert len(set(obtained)) == len(obtained)
    assert len(set(obtained)) == number_of_ids


@pytest.mark.timeout(20)
def test_get_new_request_id_threaded():
    """
    `get_new_request_id` works with multiple threads.
    """
    ids_per_run = 5
    number_of_runs = 4

    def get_many_client_id(publisher):
        return [publisher._get_new_request_id() for _ in range(ids_per_run)]

    publisher = FramePublisher()
    thread_pool = ThreadPoolExecutor(max_workers=2)
    client_id_lists = [
        thread_pool.submit(get_many_client_id, publisher).result()
        for _ in range(number_of_runs)
    ]
    obtained = list(itertools.chain(*client_id_lists))
    assert len(obtained) == ids_per_run * number_of_runs
    assert len(set(obtained)) == len(obtained)


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
@pytest.mark.parametrize("count", list(range(5)))
def test_cancellation_ends_started_stream(count):
    """
    Test that a cancelled a frame stream yields no more frames and terminates, even after consuming some frames and then
    sending more.
    """
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
