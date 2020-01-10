import time
from typing import Tuple

import pytest
from grpc import RpcError
from narupa.core.change_buffers import DictionaryChange

from narupa.core.narupa_client import NarupaClient
from narupa.core.narupa_server import NarupaServer

IMMEDIATE_REPLY_WAIT_TIME = 0.1

ACCESS_TOKEN_1 = object()
ACCESS_TOKEN_2 = object()
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
