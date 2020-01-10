import time
from typing import Tuple

import pytest
from narupa.core.change_buffers import DictionaryChange

from narupa.core.narupa_client import NarupaClient
from narupa.core.narupa_server import NarupaServer

IMMEDIATE_REPLY_WAIT_TIME = 0.1

INITIAL_STATE = {
    'hello': 100,
    'test': {'baby': 'yoda'},
}


@pytest.fixture
def client_server() -> Tuple[NarupaClient, NarupaServer]:
    with NarupaServer(address="localhost", port=0) as server:
        change = DictionaryChange(INITIAL_STATE, set())
        server.update_state(None, change)
        with NarupaClient(address="localhost", port=server.port) as client:
            yield client, server


def test_server_cant_set_non_basic_type(client_server):
    """
    Test that setting a value to a non-basic type raises a ValueError.
    """
    client, server = client_server

    class UserType:
        pass

    change = DictionaryChange({'hello': UserType()}, set())
    with pytest.raises(TypeError):
        server.update_state(None, change)


def test_client_cant_set_non_basic_type(client_server):
    """
    Test that setting a value to a non-basic type raises a ValueError.
    """
    client, server = client_server

    class UserType:
        pass

    change = DictionaryChange({'hello': UserType()}, set())
    with pytest.raises(TypeError):
        client.attempt_update_state(change)


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

    change = DictionaryChange({'hello': 'goodbye'}, {'test'})
    client.attempt_update_state(change)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    with server.lock_state() as state:
        assert state == {'hello': 'goodbye'}


def test_client_state_reflects_own_update(client_server):
    """
    Test that a client state reflects changes requested by that client.
    """
    client, server = client_server
    client.subscribe_all_state_updates(0)

    change = DictionaryChange({'hello': 'goodbye'}, {'test'})
    client.attempt_update_state(change)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    with client.lock_state() as state:
        assert state == {'hello': 'goodbye'}


def test_client_state_reflects_other_update(client_server):
    """
    Test that a client state reflects changes requested by that client.
    """
    client1, server = client_server
    client1.subscribe_all_state_updates(0)

    with NarupaClient(address="localhost", port=server.port) as client2:
        change = DictionaryChange({'hello': 'goodbye'}, {'test'})
        client2.attempt_update_state(change)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    with client1.lock_state() as state:
        assert state == {'hello': 'goodbye'}


