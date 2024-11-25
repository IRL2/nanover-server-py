import time
from typing import Tuple

import numpy
import pytest

from nanover.testing import assert_in_soon
from nanover.utilities.change_buffers import DictionaryChange
from nanover.utilities.key_lockable_map import ResourceLockedError

from nanover.core.nanover_client import NanoverClient
from nanover.core.nanover_server import NanoverServer

IMMEDIATE_REPLY_WAIT_TIME = 0.01
ARBITRARY_LOCK_DURATION = 5

INITIAL_STATE = {
    "hello": 100,
    "test": {"baby": "yoda"},
}


@pytest.fixture
def client_server() -> Tuple[NanoverClient, NanoverServer]:
    with NanoverServer(address="localhost", port=0) as server:
        change = DictionaryChange(INITIAL_STATE)
        server.update_state(None, change)
        with NanoverClient.insecure_channel(
            address="localhost", port=server.port
        ) as client:
            yield client, server


def test_server_cannot_set_non_basic_type(client_server):
    """
    Test that setting a value to a non-basic type raises a ValueError.
    """
    client, server = client_server

    class UserType:
        pass

    change = DictionaryChange({"hello": UserType()})
    with pytest.raises(TypeError):
        server.update_state(None, change)


def test_client_cannot_set_non_basic_type(client_server):
    """
    Test that setting a value to a non-basic type raises a ValueError.
    """
    client, server = client_server

    class UserType:
        pass

    change = DictionaryChange({"hello": UserType()})
    with pytest.raises(TypeError):
        client.attempt_update_state(change)


def test_client_copy_state_equals_actual_state(client_server):
    """
    Test that copying state from client gives a state equal to the client state.
    """
    client, server = client_server
    client.subscribe_all_state_updates(0)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    copy = client.copy_state()

    with client.lock_state() as state:
        assert copy == state


def test_client_copy_state_is_independent(client_server):
    """
    Test that copying state from client gives an independent state that can be
    modified without changing the client state.
    """
    client, server = client_server
    client.subscribe_all_state_updates(0)

    # wait for first update
    assert_in_soon(lambda: "test", lambda: client.copy_state())

    copy = client.copy_state()
    copy["test"]["baby"] = "shark"

    with client.lock_state() as state:
        assert copy != state


def test_server_copy_state_equals_actual_state(client_server):
    """
    Test that copying state from server gives a state equal to the server state.
    """
    client, server = client_server

    copy = server.copy_state()

    with server.lock_state() as state:
        assert copy == state


def test_server_copy_state_is_independent(client_server):
    """
    Test that copying state from server gives an independent state that can be
    modified without changing the server state.
    """
    client, server = client_server

    copy = server.copy_state()
    copy["test"]["baby"] = "shark"

    with server.lock_state() as state:
        assert copy != state


def test_server_has_initial_state(client_server):
    """
    Test that the server has the correct initial state.
    """
    client, server = client_server

    with server.lock_state() as state:
        assert state == INITIAL_STATE


def test_client_receive_initial_state(client_server):
    """
    Test that the client state matches the initial state after subscribing to
    state updates from the server.
    """
    client, server = client_server
    client.subscribe_all_state_updates(0)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    with client.lock_state() as state:
        assert state == INITIAL_STATE


def test_server_state_reflects_client_update(client_server):
    """
    Test that a server state reflects changes requested by the client.
    """
    client, server = client_server

    change = DictionaryChange({"hello": "goodbye"}, {"test"})
    client.attempt_update_state(change)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    with server.lock_state() as state:
        assert state == {"hello": "goodbye"}


def test_client_state_reflects_own_update(client_server):
    """
    Test that a client state reflects changes requested by that client.
    """
    client, server = client_server
    client.subscribe_all_state_updates(0)

    change = DictionaryChange({"hello": "goodbye"}, {"test"})
    client.attempt_update_state(change)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    with client.lock_state() as state:
        assert state == {"hello": "goodbye"}


def test_client_state_reflects_other_update(client_server):
    """
    Test that a client state reflects changes requested by that client.
    """
    client1, server = client_server
    client1.subscribe_all_state_updates(0)

    with NanoverClient.insecure_channel(
        address="localhost", port=server.port
    ) as client2:
        change = DictionaryChange({"hello": "goodbye"}, {"test"})
        client2.attempt_update_state(change)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    with client1.lock_state() as state:
        assert state == {"hello": "goodbye"}


@pytest.mark.parametrize("update_interval", (1 / 10, 1 / 30, 1 / 60))
def test_subscribe_updates_sends_initial_immediately(client_server, update_interval):
    """
    Test that subscribing updates before any have been sent will immediately
    send the initial values regardless of interval.
    """
    client, server = client_server
    client.subscribe_all_state_updates(update_interval)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    with client.lock_state() as state:
        assert state == INITIAL_STATE


@pytest.mark.parametrize("update_interval", (0.5, 0.2, 0.1))
def test_subscribe_updates_interval(client_server, update_interval):
    """
    Test that state updates are sent at the requested interval.
    """
    test_count = 8
    deadline = time.perf_counter() + update_interval * test_count + 2

    client, server = client_server
    client.subscribe_all_state_updates(update_interval)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    with client.lock_state() as state:
        assert state["hello"] == INITIAL_STATE["hello"]

    update_times = []

    for i in range(test_count):
        change = DictionaryChange({"hello": i})
        client.attempt_update_state(change)

        time_before = time.perf_counter()

        while time.perf_counter() < deadline:
            with client.lock_state() as state:
                if state["hello"] == i:
                    break
            time.sleep(0.01)
        else:
            raise Exception("Test timed out.")

        time_after = time.perf_counter()
        update_times.append(time_after - time_before)

    assert numpy.average(update_times) == pytest.approx(update_interval, abs=0.1)


def test_can_lock_unlocked(client_server):
    """
    Test that an unlocked state key can be locked.
    """
    client, server = client_server
    assert client.attempt_update_locks({"hello": ARBITRARY_LOCK_DURATION})


def test_can_lock_unlocked_nonexistent(client_server):
    """
    Test that a nonexistent unlocked state key can be locked.
    """
    client, server = client_server
    assert client.attempt_update_locks({"goodbye": ARBITRARY_LOCK_DURATION})


def test_can_lock_own_locked(client_server):
    """
    Test that an attempt to lock a state key you have already lock succeeds.
    """
    client, server = client_server
    client.attempt_update_locks({"hello": ARBITRARY_LOCK_DURATION})
    assert client.attempt_update_locks({"hello": ARBITRARY_LOCK_DURATION})


def test_can_update_own_locked(client_server):
    """
    Test that you can update state keys that you locked.
    """
    client, server = client_server
    client.attempt_update_locks({"hello": ARBITRARY_LOCK_DURATION})
    change = DictionaryChange({"hello": 999}, {"test"})
    assert client.attempt_update_state(change)


def test_can_release_own_lock(client_server):
    """
    Test that you can release your own lock.
    """
    client, server = client_server
    client.attempt_update_locks({"hello": ARBITRARY_LOCK_DURATION})
    assert client.attempt_update_locks({"hello": None})


def test_can_update_unlocked(client_server):
    """
    Test that you can update state keys that are unlocked.
    """
    client, server = client_server
    change = DictionaryChange({"hello": 999}, {"test"})
    assert client.attempt_update_state(change)


def test_cannot_remove_locked_key(client_server):
    """
    Test can't remove key that is locked by another access token.
    """
    client, server = client_server
    server.update_locks(acquire={"hello": ARBITRARY_LOCK_DURATION})
    change = DictionaryChange({}, {"hello"})
    assert not client.attempt_update_state(change)


def test_update_unlocked_repeated(client_server):
    client1, server = client_server
    change1 = DictionaryChange({"hello": 1})
    change2 = DictionaryChange({"hello": 2})
    with NanoverClient.insecure_channel(
        address="localhost", port=server.port
    ) as client2:
        assert client1.attempt_update_state(change1)
        assert client2.attempt_update_state(change2)
        assert client1.attempt_update_state(change1)
        assert client2.attempt_update_state(change2)


def test_cannot_lock_other_locked(client_server):
    """
    Test that you cannot lock a state key that is locked by someone else.
    """
    client1, server = client_server

    with NanoverClient.insecure_channel(
        address="localhost", port=server.port
    ) as client2:
        client2.attempt_update_locks({"hello": ARBITRARY_LOCK_DURATION})
        assert not client1.attempt_update_locks({"hello": ARBITRARY_LOCK_DURATION})


def test_cannot_release_other_lock(client_server):
    """
    Test that you cannot unlock a state key that is locked by someone else.
    """
    client1, server = client_server

    with NanoverClient.insecure_channel(
        address="localhost", port=server.port
    ) as client2:
        client2.attempt_update_locks({"hello": ARBITRARY_LOCK_DURATION})
        assert client1.attempt_update_locks({"hello": None})
        change = DictionaryChange({"hello": 999})
        assert not client1.attempt_update_state(change)


def test_cannot_set_other_locked(client_server):
    """
    Test that you cannot set a state key that is locked by someone else.
    """
    client1, server = client_server

    with NanoverClient.insecure_channel(
        address="localhost", port=server.port
    ) as client2:
        client2.attempt_update_locks({"hello": ARBITRARY_LOCK_DURATION})
        change = DictionaryChange({"hello": 999})
        assert not client1.attempt_update_state(change)


@pytest.mark.parametrize("lock_duration", (0.5, 1, 2))
def test_lock_durations(client_server, lock_duration):
    """
    Test that locks expire roughly after the requested duration has passed.
    """
    client, server = client_server

    access_token_1 = object()
    access_token_2 = object()

    server.update_locks(access_token_1, acquire={"hello": lock_duration})
    lock_time = time.perf_counter()
    lock_deadline = lock_time + lock_duration

    while time.perf_counter() <= lock_deadline:
        with pytest.raises(ResourceLockedError):
            server.update_locks(access_token_2, acquire={"hello": lock_duration})
        time.sleep(0.1)

    time.sleep(0.1)

    server.update_locks(access_token_2, acquire={"hello": lock_duration})
