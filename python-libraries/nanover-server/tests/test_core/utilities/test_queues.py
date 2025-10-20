"""
Unit tests for :mod:`nanover.core.request_queues`.
"""

import time
from concurrent.futures import ThreadPoolExecutor
from queue import Empty
import pytest
import itertools

from nanover.utilities.queues import LastItemQueue


class TestSingleItemQueue:
    @pytest.fixture
    def single_item_queue(self):
        return LastItemQueue()

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
            single_item_queue.get(block=True, timeout=0.5)

    def test_put_one_item(self, single_item_queue):
        item = "hello"
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
        item = "hello"
        single_item_queue.put(item)
        retrieved = single_item_queue.get(block=False)
        assert retrieved == item

    def test_many_get(self, single_item_queue):
        single_item_queue.put(0, block=False)
        single_item_queue.get(block=False)
        with pytest.raises(Empty):
            single_item_queue.get(block=False)

    @pytest.mark.timeout(20)
    @pytest.mark.parametrize("blocking", (True, False))
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
            while context["running"]:
                try:
                    # timeout only applies when blocking
                    item = queue.get(block=blocking, timeout=0.5)
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

        context = {"running": True}
        consumer_futures = [
            thread_pool.submit(consume_data, single_item_queue, context)
            for _ in range(number_of_consumers)
        ]

        # wait for the producers
        for future in producer_futures:
            future.result()

        # get the result of the consumers
        context["running"] = False
        obtained = list(
            itertools.chain(*(future.result() for future in consumer_futures))
        )

        assert len(obtained) == len(set(obtained))
        assert len(obtained) <= number_of_producers * number_of_records_per_producer
