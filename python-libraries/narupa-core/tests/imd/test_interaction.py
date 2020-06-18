import narupa.protocol.imd.imd_pb2 as imd_pb2
import numpy as np
import pytest
from narupa.imd.particle_interaction import ParticleInteraction, DEFAULT_MAX_FORCE
from narupa.utilities.protobuf_utilities import dict_to_struct
from hypothesis import strategies as st, given
from .. import *

@pytest.fixture
def interaction():
    return ParticleInteraction(player_id='test player', interaction_id='test interaction')


@pytest.fixture
def interaction_with_properties():
    return ParticleInteraction(
        player_id='test player',
        interaction_id='test interaction',
        arbitrary_property='arbitrary value',
        other_arbitrary_property='other arbitrary value',
    )


def test_player_id():
    interaction = ParticleInteraction(player_id="2", interaction_id='test interaction')
    assert interaction.player_id == "2"


def test_interaction_id():
    interaction = ParticleInteraction(player_id='test player', interaction_id="2")
    assert interaction.interaction_id == "2"


def test_get_default_position(interaction):
    assert np.allclose(interaction.position, [0, 0, 0])


def test_set_position(interaction):
    interaction.position = [1, 1, 1]
    assert np.allclose(interaction.position, [1, 1, 1])


def test_from_proto():
    interaction_grpc = imd_pb2.ParticleInteraction(
        player_id='1',
        interaction_id='0',
        position=(0, 0, 0),
    )
    interaction = ParticleInteraction.from_proto(interaction_grpc)
    assert interaction.player_id == "1"
    assert interaction.interaction_id == "0"
    assert interaction.type == "gaussian"
    assert interaction.scale == 1
    assert interaction.mass_weighted is True
    assert interaction.max_force == DEFAULT_MAX_FORCE


def test_from_proto_properties():
    struct = dict_to_struct({
        ParticleInteraction.TYPE_KEY: "harmonic",
        ParticleInteraction.SCALE_KEY: 150,
        ParticleInteraction.RESET_VELOCITIES_KEY: True,
        ParticleInteraction.MAX_FORCE_KEY: 5000,
        ParticleInteraction.MASS_WEIGHTED_KEY: False
    })
    interaction_grpc = imd_pb2.ParticleInteraction(
        player_id='1',
        interaction_id='0',
        position=(0, 0, 0),
        properties=struct,
    )
    interaction = ParticleInteraction.from_proto(interaction_grpc)
    assert interaction.player_id == "1"
    assert interaction.interaction_id == "0"
    assert interaction.type == "harmonic"
    assert interaction.scale == 150
    assert interaction.mass_weighted is False
    assert interaction.max_force == 5000
    assert interaction.reset_velocities is True


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
    proto.position[:] = (0, 0, 0)
    interaction = ParticleInteraction.from_proto(proto)
    assert interaction.mass_weighted is True


def test_set_reset_vels(interaction):
    interaction.reset_velocities = True
    assert interaction.reset_velocities is True


def test_set_mass(interaction):
    interaction.mass_weighted = False
    assert interaction.mass_weighted is False


def test_get_proto(interaction):
    proto = interaction.proto
    assert proto.player_id == interaction.player_id
    assert proto.interaction_id == interaction.interaction_id
    assert np.allclose(proto.position, [0, 0, 0])


def test_get_proto_properties(interaction_with_properties):
    proto = interaction_with_properties.proto
    assert proto.properties['arbitrary_property'] == 'arbitrary value'
    assert proto.properties['other_arbitrary_property'] == 'other arbitrary value'


@st.composite
def interactions(draw):
    position = draw(st.lists(st.floats(allow_infinity=False, max_value=MAX_FLOAT32, width=32), min_size=3, max_size=3))
    particle_ids = draw(st.lists(st.integers(min_value=0, max_value=MAX_INT32)))

    keywords = draw(st.dictionaries(st.text(), EXACT_SINGLE_VALUE_STRATEGY))

    interaction_type = draw(st.one_of(st.none(), st.text(), st.just('gaussian'), st.just('harmonic')))
    if interaction_type is not None:
        keywords['interaction_type'] = interaction_type

    scale = draw(st.one_of(st.none(), st.floats(allow_nan=False, allow_infinity=False)))
    if scale is not None:
        keywords['scale'] = scale

    reset_velocities = draw(st.one_of(st.none(), st.booleans()))
    if reset_velocities is not None:
        keywords['reset_velocities'] = reset_velocities

    mass_weighted = draw(st.one_of(st.none(), st.booleans()))
    if mass_weighted is not None:
        keywords['mass_weighted'] = mass_weighted

    max_force = draw(st.one_of(st.none(), st.floats(allow_nan=False)))
    if max_force is not None:
        keywords['max_force'] = max_force

    player_id = draw(st.text())
    interaction_id = draw(st.text())

    return ParticleInteraction(player_id, interaction_id, position, particle_ids, **keywords)


@given(interactions())
def test_serialize_then_deserialize(interaction):
    proto = interaction.proto
    new_interaction = ParticleInteraction.from_proto(proto)
    assert new_interaction == interaction
