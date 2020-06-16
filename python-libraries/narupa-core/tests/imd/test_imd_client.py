import time
import grpc
import pytest

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


def test_subscribe_interactions(imd_server_client, interaction):
    """
    Test that IMD interactions can be subscribed.
    """
    imd_server, imd_client = imd_server_client
    imd_client.subscribe_interactions()


def test_subscribe_own_interaction(imd_server_client, interaction):
    """
    Test that after subscribing interactions, we receive our own published
    interaction.
    """
    imd_server, imd_client = imd_server_client
    imd_client.subscribe_interactions()
    local_interaction_id = imd_client.start_interaction()
    imd_client.update_interaction(local_interaction_id, interaction)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME * 5)
    assert interaction.interaction_id in imd_client.interactions


def test_subscribe_own_interaction_removed(imd_server_client, interaction):
    """
    Test that after subscribing interactions, we receive our own published
    interaction and after removal it is removed.
    """
    imd_server, imd_client = imd_server_client
    imd_client.subscribe_interactions()
    local_interaction_id = imd_client.start_interaction()
    imd_client.update_interaction(local_interaction_id, interaction)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    assert interaction.interaction_id in imd_client.interactions

    imd_client.stop_interaction(local_interaction_id)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)
    assert interaction.interaction_id not in imd_client.interactions
