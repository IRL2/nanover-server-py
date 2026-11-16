"""
Unit tests of the IMD service, without any connections.
"""

from unittest.mock import Mock

from nanover.imd.imd_state import (
    VELOCITY_RESET_KEY,
    ImdStateWrapper,
    interaction_to_dict,
)
from nanover.imd.particle_interaction import ParticleInteraction
from nanover.utilities.change_buffers import DictionaryChange
from nanover.utilities.state_dictionary import StateDictionary


def test_add_duplicate_interaction_id():
    imd_state = ImdStateWrapper(StateDictionary())
    imd_state.insert_interaction("interaction.test", ParticleInteraction())
    imd_state.insert_interaction("interaction.test", ParticleInteraction())
    assert len(imd_state.active_interactions) == 1


def test_multiple_keys():
    imd_state = ImdStateWrapper(StateDictionary())
    imd_state.insert_interaction("interaction.test1", ParticleInteraction())
    imd_state.insert_interaction("interaction.test2", ParticleInteraction())
    assert len(imd_state.active_interactions) == 2


def test_velocity_reset_enabled():
    state = StateDictionary()
    imd_state = ImdStateWrapper(state)
    imd_state.velocity_reset_available = True
    assert imd_state.velocity_reset_available
    assert state.copy_content()[VELOCITY_RESET_KEY]


def test_interaction_started_event_insert():
    imd_state = ImdStateWrapper(StateDictionary())

    mock = Mock()
    imd_state.interaction_started.add_callback(mock)

    id = "interaction.test1"
    interaction = ParticleInteraction()

    imd_state.insert_interaction(id, interaction)
    imd_state.insert_interaction(id, interaction)
    imd_state.insert_interaction(id, interaction)

    mock.assert_called_once_with(key=id, interaction=interaction)


def test_interaction_started_event_raw():
    state = StateDictionary()
    imd_state = ImdStateWrapper(state)

    mock = Mock()
    imd_state.interaction_started.add_callback(mock)

    id = "interaction.test1"
    interaction = ParticleInteraction()

    state.update_state(
        None, DictionaryChange(updates={id: interaction_to_dict(interaction)})
    )

    mock.assert_called_once_with(key=id, interaction=interaction)


def test_interaction_stopped_event():
    imd_state = ImdStateWrapper(StateDictionary())

    mock = Mock()
    imd_state.interaction_stopped.add_callback(mock)

    id = "interaction.test1"
    interaction = ParticleInteraction()

    imd_state.insert_interaction(id, interaction)
    imd_state.remove_interaction(id)
    imd_state.remove_interaction(id)
    imd_state.remove_interaction(id)

    mock.assert_called_once_with(key=id, interaction=interaction)


def test_interaction_stopped_event_raw():
    state = StateDictionary()
    imd_state = ImdStateWrapper(state)

    mock = Mock()
    imd_state.interaction_stopped.add_callback(mock)

    id = "interaction.test1"
    interaction = ParticleInteraction()

    imd_state.insert_interaction(id, interaction)
    state.update_state(None, DictionaryChange(removals={id}))

    mock.assert_called_once_with(key=id, interaction=interaction)
