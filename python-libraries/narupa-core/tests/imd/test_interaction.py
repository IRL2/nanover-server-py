import numpy as np
import pytest

from narupa.imd.interaction import Interaction


@pytest.fixture
def interaction():
    return Interaction()


def test_player_id():
    interaction = Interaction("2")
    assert interaction.player_id == "2"

def test_interaction_id():
    interaction = Interaction(interaction_id = "2")
    assert interaction.interaction_id == "2"

def test_get_default_position(interaction):
    assert np.allclose(interaction.position, [0, 0, 0])


def test_set_position(interaction):
    interaction.position = [1, 1, 1]
    assert np.allclose(interaction.position, [1, 1, 1])


def test_set_invalid_position(interaction):
    with pytest.raises(ValueError):
        interaction.position = [0, 0]


def test_get_default_particles(interaction):
    assert len(interaction.particles) == 0


def test_set_particles(interaction):
    interaction.particles = [0, 1, 2, 3, 4]
    assert np.allclose(interaction.particles, [0, 1, 2, 3, 4])


def test_set_particle_unique(interaction):
    interaction.particles = [0, 0, 0, 1, 2, 3, 4]
    assert np.allclose(interaction.particles, [0, 1, 2, 3, 4])

def test_set_property_number(interaction):
    interaction.properties['property'] = 2.0
    assert interaction.properties['property'] == pytest.approx(2.0)

def test_set_property_str(interaction):
    interaction.properties['property'] = 'value'
    assert interaction.properties['property'] == 'value'

def test_set_property_list(interaction):
    interaction.properties['property'] = [5, 4, 3, 2, 1]
    assert np.allclose(interaction.properties['property'], [5,4,3,2,1])
