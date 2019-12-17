import time
from concurrent import futures
from queue import Queue

import grpc
import pytest
from narupa.imd.imd_client import queue_generator

from .test_imd_server import imd_server_client, imd_server, interaction

IMMEDIATE_REPLY_WAIT_TIME = 0.01


def test_start_interaction(imd_server_client):
    """
    Test that you can start an interaction.
    """
    imd_server, imd_client = imd_server_client
    assert imd_client.start_interaction() is not None


def test_start_interaction_twice(imd_server_client):
    """
    Test that you can start two independent interactions.
    """
    imd_server, imd_client = imd_server_client
    id1 = imd_client.start_interaction()
    id2 = imd_client.start_interaction()
    assert id1 != id2


def test_update_interaction(imd_server_client, interaction):
    imd_server, imd_client = imd_server_client
    imd_client.start_interaction()
    interaction_id = imd_client.start_interaction()

    imd_client.update_interaction(interaction_id, interaction)


def test_update_unknown_interaction(imd_server_client, interaction):
    imd_server, imd_client = imd_server_client
    imd_client.start_interaction()
    interaction_id = imd_client.start_interaction()

    with pytest.raises(KeyError):
        imd_client.update_interaction(interaction_id + "nonsense", interaction)


def test_delete_interaction(imd_server_client, interaction):
    imd_server, imd_client = imd_server_client
    interaction_id = imd_client.start_interaction()
    imd_client.stop_interaction(interaction_id)
    assert len(imd_client._local_interactions_states) == 0


def test_delete_unknown_interaction(imd_server_client, interaction):
    imd_server, imd_client = imd_server_client
    interaction_id = imd_client.start_interaction()
    with pytest.raises(KeyError):
        imd_client.stop_interaction(interaction_id + "nonsense")


def test_delete_deleted_interaction(imd_server_client, interaction):
    imd_server, imd_client = imd_server_client
    interaction_id = imd_client.start_interaction()
    imd_client.stop_interaction(interaction_id)
    with pytest.raises(KeyError):
        imd_client.stop_interaction(interaction_id)


def test_update_deleted_interaction(imd_server_client, interaction):
    imd_server, imd_client = imd_server_client
    interaction_id = imd_client.start_interaction()
    imd_client.stop_interaction(interaction_id)
    with pytest.raises(KeyError):
        imd_client.update_interaction(interaction_id, interaction)


def test_stop_all_interactions(imd_server_client, interaction):
    imd_server, imd_client = imd_server_client
    imd_client.start_interaction()
    imd_client.start_interaction()
    imd_client.stop_all_interactions()
    assert len(imd_client._local_interactions_states) == 0


def test_bad_interaction_type(imd_server_client):
    imd_server, imd_client = imd_server_client
    interaction_id = imd_client.start_interaction()
    imd_client.update_interaction(interaction_id, "something_stupid")
    with pytest.raises(grpc.RpcError):
        imd_client.stop_interaction(interaction_id)


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


def test_subscribe_interactions(imd_server_client, interaction):
    """
    Test that IMD interactions can be subscribed.
    """
    imd_server, imd_client = imd_server_client
    imd_client.subscribe_interactions(interval=0)


def test_subscribe_own_interaction(imd_server_client, interaction):
    """
    Test that after subscribing interactions, we receive our own published
    interaction.
    """
    imd_server, imd_client = imd_server_client
    imd_client.subscribe_interactions(interval=0)
    local_interaction_id = imd_client.start_interaction()
    imd_client.update_interaction(local_interaction_id, interaction)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME * 5)
    assert local_interaction_id in imd_client.interactions


def test_subscribe_own_interaction_removed(imd_server_client, interaction):
    """
    Test that after subscribing interactions, we receive our own published
    interaction and after removal it is removed.
    """
    imd_server, imd_client = imd_server_client
    imd_client.subscribe_interactions(interval=0)
    local_interaction_id = imd_client.start_interaction()
    imd_client.update_interaction(local_interaction_id, interaction)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    assert local_interaction_id in imd_client.interactions

    imd_client.stop_interaction(local_interaction_id)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    assert local_interaction_id not in imd_client.interactions


@pytest.mark.parametrize('update_interval', (1/10, 1/30, 1/60))
def test_subscribe_interactions_sends_initial_immediately(
        imd_server_client,
        update_interval,
        interaction,
):
    """
    Test that subscribing interactions before any have been sent will
    immediately send the first update regardless of interval.
    """
    imd_server, imd_client = imd_server_client
    imd_client.subscribe_interactions(interval=update_interval)
    local_interaction_id = imd_client.start_interaction()
    imd_client.update_interaction(local_interaction_id, interaction)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    assert local_interaction_id in imd_client.interactions


@pytest.mark.parametrize('update_interval', (.5, .2, .1))
def test_subscribe_interactions_interval(
        imd_server_client,
        update_interval,
        interaction,
):
    """
    Test that interaction updates are sent at the requested interval.
    """
    imd_server, imd_client = imd_server_client
    imd_client.subscribe_interactions(interval=update_interval)
    local_interaction_id = imd_client.start_interaction()

    interaction.position = [1, 0, 0]
    imd_client.update_interaction(local_interaction_id, interaction)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    interaction.position = [2, 0, 0]
    imd_client.update_interaction(local_interaction_id, interaction)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    assert imd_client.interactions[local_interaction_id].position[0] == 1

    imd_client.update_interaction(local_interaction_id, interaction)
    time.sleep(update_interval)
    assert imd_client.interactions[local_interaction_id].position[0] == 2
