"""
Unit tests of the IMD service, without any connections.
"""

from narupa.imd.imd_state import ImdStateWrapper, VELOCITY_RESET_KEY
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.state.state_dictionary import StateDictionary


def test_add_duplicate_interaction_id():
    imd_state = ImdStateWrapper(StateDictionary())
    imd_state.insert_interaction('interaction.test', ParticleInteraction())
    imd_state.insert_interaction('interaction.test', ParticleInteraction())
    assert len(imd_state.active_interactions) == 1


def test_multiple_keys():
    imd_state = ImdStateWrapper(StateDictionary())
    imd_state.insert_interaction('interaction.test1', ParticleInteraction())
    imd_state.insert_interaction('interaction.test2', ParticleInteraction())
    assert len(imd_state.active_interactions) == 2


def test_velocity_reset_enabled():
    state = StateDictionary()
    imd_state = ImdStateWrapper(state)
    imd_state.velocity_reset_available = True
    assert imd_state.velocity_reset_available
    assert state.copy_content()[VELOCITY_RESET_KEY]
