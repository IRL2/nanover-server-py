"""
Unit tests of the IMD service, without any connections.
"""

import pytest

from narupa.imd.imd_service import ImdService
from narupa.imd.particle_interaction import ParticleInteraction


@pytest.fixture
def interaction():
    return ParticleInteraction()


def test_get_key(interaction):
    key = ImdService.get_key(interaction)
    assert key == ("1", "0")


def test_add_same_key(interaction):
    service = ImdService()
    key = ImdService.get_key(interaction)
    service.interactions[key] = interaction
    interaction = ParticleInteraction()
    key = ImdService.get_key(interaction)
    service.interactions[key] = interaction
    assert len(service.interactions) == 1


def test_multiple_keys(interaction):
    service = ImdService()
    key = ImdService.get_key(interaction)
    service.interactions[key] = interaction
    interaction = ParticleInteraction()
    key = ("2", "0")
    service.interactions[key] = interaction
    assert len(service.interactions) == 2
