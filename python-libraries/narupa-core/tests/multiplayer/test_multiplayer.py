# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Integration tests of the multiplayer server with the reference multiplayer client.
"""
import pytest

from narupa.multiplayer.multiplayer_client import MultiplayerClient
from narupa.multiplayer.multiplayer_server import MultiplayerServer


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


def test_join_multiplayer(server_client_pair):
    """
    Test that it's possible to join multiplayer and receive a player id.
    """
    server, client = server_client_pair
    player_id = client.create_player_id()
    assert player_id is not None


def test_join_multiplayer_twice_same_id(server_client_pair):
    """
    Test that joining multiplayer again gives you your existing player id.
    """
    server, client = server_client_pair
    first_id = client.create_player_id()
    second_id = client.create_player_id()
    assert first_id == second_id
