"""
Unit tests of the IMD service, without any connections.
"""

import pytest

from narupa.imd.imd_service import ImdService
from narupa.imd.particle_interaction import ParticleInteraction


@pytest.fixture
def interaction():
    return ParticleInteraction(player_id='test player', interaction_id='test interaction')


def test_add_duplicate_interaction_id(interaction):
    service = ImdService()
    service.insert_interaction(interaction)
    interaction = ParticleInteraction(
        player_id=interaction.player_id,
        interaction_id=interaction.interaction_id,
    )
    service.insert_interaction(interaction)
    assert len(service.active_interactions) == 1


def test_multiple_keys(interaction):
    interaction2 = ParticleInteraction(player_id="T", interaction_id="T.0")

    service = ImdService()
    service.insert_interaction(interaction)
    service.insert_interaction(interaction2)
    assert len(service.active_interactions) == 2
