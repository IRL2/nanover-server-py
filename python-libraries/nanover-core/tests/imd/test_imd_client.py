import itertools
import time
from typing import Callable

import pytest
from nanover.imd.particle_interaction import ParticleInteraction
from nanover.imd.imd_state import (
    INTERACTION_PREFIX,
    interaction_to_dict,
    dict_to_interaction,
)
from nanover.utilities.change_buffers import DictionaryChange

from .test_imd_server import imd_server_client, imd_server, interaction

IMMEDIATE_REPLY_WAIT_TIME = 0.01


def assert_equal_soon(a: Callable, b: Callable, interval=0.1, timeout=1.0):
    __tracebackhide__ = True
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline and a() != b():
        time.sleep(interval)
    assert a() == b()


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


def test_update_interaction(imd_server_client):
    imd_server, imd_client = imd_server_client
    imd_client.start_interaction()
    interaction_id = imd_client.start_interaction()
    interaction = ParticleInteraction()

    imd_client.update_interaction(interaction_id, interaction)


def test_update_unknown_interaction(imd_server_client):
    imd_server, imd_client = imd_server_client
    imd_client.start_interaction()
    interaction_id = imd_client.start_interaction()
    interaction = ParticleInteraction()

    with pytest.raises(KeyError):
        imd_client.update_interaction(interaction_id + "nonsense", interaction)


def test_delete_interaction(imd_server_client):
    imd_server, imd_client = imd_server_client
    interaction_id = imd_client.start_interaction()
    imd_client.stop_interaction(interaction_id)
    assert len(imd_client._local_interaction_ids) == 0


def test_delete_unknown_interaction(imd_server_client):
    imd_server, imd_client = imd_server_client
    interaction_id = imd_client.start_interaction()
    with pytest.raises(KeyError):
        imd_client.stop_interaction(interaction_id + "nonsense")


def test_delete_deleted_interaction(imd_server_client):
    imd_server, imd_client = imd_server_client
    interaction_id = imd_client.start_interaction()
    imd_client.stop_interaction(interaction_id)
    with pytest.raises(KeyError):
        imd_client.stop_interaction(interaction_id)


def test_update_deleted_interaction(imd_server_client):
    imd_server, imd_client = imd_server_client
    interaction_id = imd_client.start_interaction()
    interaction = ParticleInteraction()
    imd_client.stop_interaction(interaction_id)
    with pytest.raises(KeyError):
        imd_client.update_interaction(interaction_id, interaction)


def test_stop_all_interactions(imd_server_client):
    imd_server, imd_client = imd_server_client
    imd_client.start_interaction()
    imd_client.start_interaction()
    imd_client.stop_all_interactions()
    assert len(imd_client._local_interaction_ids) == 0


def test_bad_interaction_type(imd_server_client):
    imd_server, imd_client = imd_server_client
    interaction_id = imd_client.start_interaction()
    with pytest.raises(TypeError):
        imd_client.update_interaction(interaction_id, "something_stupid")


def test_subscribe_interactions(imd_server_client):
    """
    Test that IMD interactions can be subscribed.
    """
    imd_server, imd_client = imd_server_client
    imd_client.subscribe_all_state_updates(interval=0)


def test_subscribe_own_interaction(imd_server_client):
    """
    Test that after subscribing interactions, we receive our own published
    interaction.
    """
    imd_server, imd_client = imd_server_client
    imd_client.subscribe_all_state_updates(interval=0)

    interaction_id = imd_client.start_interaction()
    interaction = ParticleInteraction()
    imd_client.update_interaction(interaction_id, interaction)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME * 5)
    assert interaction_id in imd_client.interactions


def test_subscribe_own_interaction_removed(imd_server_client):
    """
    Test that after subscribing interactions, we receive our own published
    interaction and after removal it is removed.
    """
    imd_server, imd_client = imd_server_client
    imd_client.subscribe_all_state_updates(interval=0)

    interaction_id = imd_client.start_interaction()
    interaction = ParticleInteraction()

    imd_client.update_interaction(interaction_id, interaction)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    assert interaction_id in imd_client.interactions
    assert imd_client.stop_interaction(interaction_id)

    time.sleep(IMMEDIATE_REPLY_WAIT_TIME * 5)
    assert interaction_id not in imd_client.interactions


def test_interactions_property(imd_server_client):
    """
    Test that the `interactions` property returns what is expected.
    """
    imd_server, imd_client = imd_server_client

    # The IMD server does some validation and conversion of the interactions.
    # We want to add invalid interactions in the shared state to emulate a
    # possible third party server or a bug in the server. So we need to remove
    # the conversion step.
    imd_server._state_service.state_dictionary.content_updated._callbacks = []

    # There are the interactions we should get when calling the property.
    real_interactions = {
        f"{INTERACTION_PREFIX}.first_id": ParticleInteraction(),
        f"{INTERACTION_PREFIX}.second_id": ParticleInteraction(),
    }
    # We should not get these interactions because the ID is not one of an
    # interaction.
    interactions_with_incompatible_id = {
        "not.a.valid.interaction.id": ParticleInteraction(),
    }
    interaction_updates = DictionaryChange(
        updates={
            key: interaction_to_dict(interaction)
            for key, interaction in itertools.chain(
                real_interactions.items(),
                interactions_with_incompatible_id.items(),
            )
        }
    )
    imd_server.update_state(None, interaction_updates)

    # These interactions have an ID matching an interaction, but the content does not match.
    incorrect_interactions = {
        f"{INTERACTION_PREFIX}.third_id": "not.an.interaction",
    }
    interaction_updates = DictionaryChange(updates=incorrect_interactions)
    imd_server.update_state(None, interaction_updates)

    imd_client.subscribe_all_state_updates(interval=0)
    assert_equal_soon(
        lambda: imd_client.interactions.keys(), lambda: real_interactions.keys()
    )
