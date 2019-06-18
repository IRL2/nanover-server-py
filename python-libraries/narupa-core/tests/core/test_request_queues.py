"""
Unit tests for :mod:`narupa.core.request_queues`.
"""

import time
from concurrent.futures import ThreadPoolExecutor
import pytest

from narupa.core import request_queues


def test_one_queue_serial():
    many_queues = request_queues.DictOfQueues()
    assert not many_queues.queues
    with many_queues.one_queue(0) as queue:
        assert list(many_queues.queues.keys()) == [0]
    assert not many_queues.queues


@pytest.mark.timeout(20)
def test_one_queue_threaded():
    max_workers = 2

    def use_queues(queue_dict, thread_index, number_of_queues):
        for queue_index in range(number_of_queues):
            time.sleep(0.01)
            request_id = (thread_index, queue_index)
            with queue_dict.one_queue(request_id) as queue:
                queue.put(0)
                with queue_dict.lock:
                    assert len(queue_dict.queues) <= max_workers

    many_dict = request_queues.DictOfQueues()
    thread_pool = ThreadPoolExecutor(max_workers=max_workers)

    for thread_index in range(max_workers):
        thread_pool.submit(use_queues, many_dict, thread_index, 10)


def test_iter_queues_serial():
    many_queues = request_queues.DictOfQueues()
    number_of_queues = 10

    # This is not the recommended way to populate the dictionary of queues: it
    # is not thread safe!
    for request_id in range(number_of_queues):
        # The values should be queues, but strings are easier to compare. The
        # values are strings and not ints to differentiate them from the keys.
        many_queues.queues[request_id] = str(request_id)

    obtained = [value for value in many_queues.iter_queues()]
    expected = [str(request_id) for request_id in range(number_of_queues)]

    assert obtained == expected


@pytest.mark.timeout(20)
def test_iter_queues_threaded():
    # There is no assertion in this test. We are making sure that their is no
    # exception raised.
    def register_and_unregister_queues(queue_dict, thread_index, number_of_queues):
        for queue_index in range(number_of_queues):
            request_id = (thread_index, queue_index)
            with queue_dict.one_queue(request_id):
                time.sleep(0.01)

    def iterate_over_queues(queue_dict):
        for _ in range(10):
            for queue in queue_dict.iter_queues():
                queue.put(0)
            time.sleep(0.01)

    many_queues = request_queues.DictOfQueues()
    max_workers = 2
    thread_pool = ThreadPoolExecutor(max_workers=max_workers)
    for thread_index in range(max_workers):
        thread_pool.submit(
            register_and_unregister_queues,
            many_queues, thread_index, 10,
        )
        thread_pool.submit(iterate_over_queues, many_queues)

