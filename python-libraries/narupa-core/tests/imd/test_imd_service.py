"""
Unit tests of the IMD service, without any connections.
"""

import pytest

from narupa.imd.imd_service import ImdService
from narupa.imd.particle_interaction import ParticleInteraction


@pytest.fixture
def interaction():
    return ParticleInteraction()


def test_add_duplicate_interaction_id(interaction):
    service = ImdService()
    service.insert_interaction(interaction)
    interaction = ParticleInteraction()
    service.insert_interaction(interaction)
    assert len(service.active_interactions) == 1


def test_multiple_keys(interaction):
    interaction2 = ParticleInteraction(player_id="T", interaction_id="T.0")

    service = ImdService()
    service.insert_interaction(interaction)
    service.insert_interaction(interaction2)
    assert len(service.active_interactions) == 2
