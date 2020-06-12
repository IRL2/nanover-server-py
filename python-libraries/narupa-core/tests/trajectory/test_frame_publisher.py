"""
Unit tests for `narupa.trajectory.frame_publisher
"""

from concurrent.futures import ThreadPoolExecutor
import itertools
import pytest

from narupa.trajectory import FramePublisher
from narupa.trajectory.frame_data import FrameData, SERVER_TIMESTAMP
from narupa.protocol.trajectory import FrameData as RawFrameData


def test_user_queue():
    """
    The `_user_queue` works as expected for a context manager.
    """
    publisher = FramePublisher()
    assert not publisher.frame_queues.queues
    with publisher.frame_queues.one_queue(0):
        assert list(publisher.frame_queues.queues.keys()) == [0]
    assert not publisher.frame_queues.queues


def test_send_wrapped_frame_data():
    publisher = FramePublisher()
    frame = FrameData()
    publisher.send_frame(3, frame)
    assert publisher.last_frame_index == 3
    assert SERVER_TIMESTAMP in publisher.last_frame.values


def test_send_raw_frame_data():
    publisher = FramePublisher()
    frame = RawFrameData()
    publisher.send_frame(5, frame)
    assert publisher.last_frame_index == 5
    assert SERVER_TIMESTAMP in publisher.last_frame.values


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
