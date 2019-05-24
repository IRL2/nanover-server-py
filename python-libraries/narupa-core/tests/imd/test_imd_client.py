from concurrent import futures
from queue import Queue

import pytest
from narupa.imd.imd_client import queue_generator
from .test_imd_server import imd_client, interaction


def test_start_interaction(imd_client):
    interaction_id = imd_client.start_interaction()
    assert interaction_id == 0


def test_start_interaction_twice(imd_client):
    imd_client.start_interaction()
    interaction_id = imd_client.start_interaction()
    assert interaction_id == 1


def test_update_interaction(imd_client, interaction):
    imd_client.start_interaction()
    interaction_id = imd_client.start_interaction()

    imd_client.update_interaction(interaction_id, interaction)


def test_update_unkown_interaction(imd_client, interaction):
    imd_client.start_interaction()
    interaction_id = imd_client.start_interaction()

    with pytest.raises(KeyError):
        imd_client.update_interaction(interaction_id + 1, interaction)


def test_delete_interaction(imd_client, interaction):
    interaction_id = imd_client.start_interaction()
    imd_client.stop_interaction(interaction_id)
    assert len(imd_client._active_interactions) == 0


def test_delete_unknown_interaction(imd_client, interaction):
    interaction_id = imd_client.start_interaction()
    with pytest.raises(KeyError):
        imd_client.stop_interaction(interaction_id + 1)


def test_delete_deleted_interaction(imd_client, interaction):
    interaction_id = imd_client.start_interaction()
    imd_client.stop_interaction(interaction_id)
    with pytest.raises(KeyError):
        imd_client.stop_interaction(interaction_id)


def test_update_deleted_interaction(imd_client, interaction):
    interaction_id = imd_client.start_interaction()
    imd_client.stop_interaction(interaction_id)
    with pytest.raises(KeyError):
        imd_client.update_interaction(interaction_id, interaction)


def test_queue_generator():
    queue = Queue()
    sentinel = object()
    items = [x for x in range(10)]
    for i in items:
        queue.put(i)
    queue.put(sentinel)
    result = [x for x in queue_generator(queue, sentinel)]
    assert result == items

def to_list(generator):
    return [x for x in generator]

def test_queue_generator_threaded():
    """
    tests that running the queue generator in a thread produces the expected results.
    """
    queue = Queue()
    sentinel = object()
    threads = futures.ThreadPoolExecutor(max_workers=10)
    generator = queue_generator(queue, sentinel)
    future = threads.submit(to_list, generator)

    # submit items to the queue, which will be proceed in the other thread.
    items = [x for x in range(10)]
    for i in items:
        queue.put(i)

    queue.put(sentinel)

    result = future.result(timeout=0.01)
    assert result == items
