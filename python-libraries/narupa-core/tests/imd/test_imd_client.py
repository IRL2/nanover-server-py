import pytest

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
