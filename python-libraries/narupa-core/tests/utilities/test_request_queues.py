"""
Unit tests for :mod:`narupa.core.request_queues`.
"""

import time
from concurrent.futures import ThreadPoolExecutor
from queue import Queue, Empty
import pytest
import itertools

from narupa.utilities import request_queues


def test_one_queue_serial():
    many_queues = request_queues.DictOfQueues()
    assert not many_queues.queues
    with many_queues.one_queue(0) as queue:
        assert list(many_queues.queues.keys()) == [0]
    assert not many_queues.queues


@pytest.mark.parametrize('queue_type', (Queue, request_queues.SingleItemQueue))
def test_one_queue_type(queue_type):
    many_queues = request_queues.DictOfQueues()
    with many_queues.one_queue(0, queue_class=queue_type) as queue:
        assert isinstance(queue, queue_type)


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


class TestSingleItemQueue:
    @pytest.fixture
    def single_item_queue(self):
        return request_queues.SingleItemQueue()

    def test_queue_none(self, single_item_queue):
        single_item_queue.put(None)
        assert single_item_queue.get() is None

    @pytest.mark.timeout(3)
    def test_blocking_get_with_content(self, single_item_queue):
        single_item_queue.put(0)
        assert single_item_queue.get(block=True) == 0

    @pytest.mark.timeout(3)
    def test_blocking_get_timeout(self, single_item_queue):
        with pytest.raises(Empty):
            single_item_queue.get(block=True, timeout=.5)

    def test_put_one_item(self, single_item_queue):
        item = 'hello'
        single_item_queue.put(item)
        assert single_item_queue._item is item

    def test_put_many_item_serial(self, single_item_queue):
        for item in range(5):
            single_item_queue.put(item)
        assert single_item_queue._item is item

    @pytest.mark.timeout(20)
    def test_put_many_item_threaded(self, single_item_queue):

        def put_values(thread_id, queue):
            for i in range(10):
                time.sleep(0.01)
                item = (thread_id, i)
                queue.put(item)

        max_workers = 2
        thread_pool = ThreadPoolExecutor(max_workers=max_workers)
        futures = [
            thread_pool.submit(put_values, thread_id, single_item_queue)
            for thread_id in range(10)
        ]

        # wait for the threads to be done
        for future in futures:
            future.result()

        assert single_item_queue._item is not None

    def test_get_initial(self, single_item_queue):
        with pytest.raises(Empty):
            single_item_queue.get(block=False)

    def test_get_item(self, single_item_queue):
        item = 'hello'
        single_item_queue.put(item)
        retrieved = single_item_queue.get(block=False)
        assert retrieved == item

    def test_many_get(self, single_item_queue):
        single_item_queue.put(0, block=False)
        single_item_queue.get(block=False)
        with pytest.raises(Empty):
            single_item_queue.get(block=False)

    @pytest.mark.timeout(20)
    @pytest.mark.parametrize('blocking', (True, False))
    def test_put_and_get_threaded(self, single_item_queue, blocking):
        """
        Add and get data from the SingleItemQueue from multiple threads.
        """

        def produce_data(thread_id, queue, number_of_records):
            for i in range(number_of_records):
                time.sleep(0.01)
                item = (thread_id, i)
                queue.put(item)

        def consume_data(queue, context):
            obtained_data = []
            while context['running']:
                try:
                    # timeout only applies when blocking
                    item = queue.get(block=blocking, timeout=.5)
                except Empty:
                    pass
                else:
                    obtained_data.append(item)
            return obtained_data

        number_of_producers = 10
        number_of_consumers = 5
        number_of_records_per_producer = 10
        max_workers = 2
        thread_pool = ThreadPoolExecutor(max_workers=max_workers)
        producer_futures = [
            thread_pool.submit(
                produce_data,
                thread_id,
                single_item_queue,
                number_of_records_per_producer,
            )
            for thread_id in range(number_of_producers)
        ]

        context = {'running': True}
        consumer_futures = [
            thread_pool.submit(consume_data, single_item_queue, context)
            for _ in range(number_of_consumers)
        ]

        # wait for the producers
        for future in producer_futures:
            future.result()

        # get the result of the consumers
        context['running'] = False
        obtained = list(itertools.chain(*(
            future.result() for future in consumer_futures
        )))

        assert len(obtained) == len(set(obtained))
        assert len(obtained) <= number_of_producers * number_of_records_per_producer
