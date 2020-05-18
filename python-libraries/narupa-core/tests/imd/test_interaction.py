import narupa.protocol.imd.imd_pb2 as imd_pb2
import numpy as np
import pytest
from narupa.imd.particle_interaction import ParticleInteraction, DEFAULT_MAX_FORCE


@pytest.fixture
def interaction():
    return ParticleInteraction(interaction_id='invalid id')


def test_interaction_id():
    interaction = ParticleInteraction(interaction_id="2")
    assert interaction.interaction_id == "2"


def test_get_default_position(interaction):
    assert np.allclose(interaction.position, [0, 0, 0])


def test_set_position(interaction):
    interaction.position = [1, 1, 1]
    assert np.allclose(interaction.position, [1, 1, 1])


def test_from_proto():
    interaction_grpc = imd_pb2.ParticleInteraction(interaction_id='0')
    interaction = ParticleInteraction.from_proto(interaction_grpc)
    assert interaction.interaction_id == "0"
    assert interaction.type == "gaussian"
    assert interaction.scale == 1
    assert interaction.mass_weighted is True
    assert interaction.max_force == DEFAULT_MAX_FORCE


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
    assert np.allclose(interaction.properties['property'], [5, 4, 3, 2, 1])


def test_get_type(interaction):
    assert interaction.type == "gaussian"


def test_set_type(interaction):
    interaction.type = "harmonic"
    assert interaction.type == "harmonic"


def test_get_scale(interaction):
    assert interaction.scale == 1


def test_set_scale(interaction):
    interaction.scale = 2
    assert interaction.scale == 2


def test_get_mass(interaction):
    assert interaction.mass_weighted is True


def test_get_mass_unset():
    proto = imd_pb2.ParticleInteraction()
    interaction = ParticleInteraction.from_proto(proto)
    assert interaction.mass_weighted is True


def test_set_reset_vels(interaction):
    interaction.reset_velocities = True
    assert interaction.properties['reset_velocities'] is True


def test_set_mass(interaction):
    interaction.mass_weighted = False
    assert interaction.properties['mass_weighted'] is False


def test_get_proto(interaction):
    proto = interaction.proto
    assert proto.interaction_id == interaction.interaction_id
    assert np.allclose(proto.position, [0, 0, 0])
