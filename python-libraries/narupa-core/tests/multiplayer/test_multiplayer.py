# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Integration tests of the multiplayer server with the reference multiplayer client.
"""
import pytest

from narupa.multiplayer.multiplayer_client import MultiplayerClient
from narupa.multiplayer.multiplayer_server import MultiplayerServer

CONNECT_WAIT_TIME = 0.01
IMMEDIATE_REPLY_WAIT_TIME = 0.01


@pytest.fixture
def server_client_pair():
    """
    Provides a multiplayer server hosting on an available port on localhost,
    and a multiplayer client connected to it.
    """
    server = MultiplayerServer(address='localhost', port=0)
    client = MultiplayerClient.insecure_channel(port=server.port)

    with client, server:
        yield server, client


@pytest.fixture
def scene():
    """
    Provides scene test data.
    """
    return {
        'position': {"x": 1., "y": 1., "z": 1.},
        'rotation': {"x": 0., "y": 0., "z": 0., "w": 1.},
        'scale': 1.,
    }


def test_join_multiplayer(server_client_pair):
    """
    Test that it's possible to join multiplayer and receive a player id.
    """
    server, client = server_client_pair
    player_id = client.join_multiplayer("user")
    assert player_id is not None


def test_join_multiplayer_twice_same_id(server_client_pair):
    """
    Test that joining multiplayer again gives you your existing player id.
    """
    server, client = server_client_pair
    first_id = client.join_multiplayer("user")
    second_id = client.join_multiplayer("user")
    assert first_id == second_id
