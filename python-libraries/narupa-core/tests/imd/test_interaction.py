import pytest
from narupa.imd.particle_interaction import ParticleInteraction

from .. import *

@pytest.fixture
def interaction():
    return ParticleInteraction()


@pytest.fixture
def interaction_with_properties():
    return ParticleInteraction(
        arbitrary_property='arbitrary value',
        other_arbitrary_property='other arbitrary value',
    )


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
    assert np.allclose(interaction.properties['property'], [5, 4, 3, 2, 1])


def test_get_type(interaction):
    assert interaction.interaction_type == "gaussian"


def test_set_type(interaction):
    interaction.interaction_type = "harmonic"
    assert interaction.interaction_type == "harmonic"


def test_get_scale(interaction):
    assert interaction.scale == 1


def test_set_scale(interaction):
    interaction.scale = 2
    assert interaction.scale == 2


def test_get_mass(interaction):
    assert interaction.mass_weighted is True


def test_set_reset_vels(interaction):
    interaction.reset_velocities = True
    assert interaction.reset_velocities is True


def test_set_mass(interaction):
    interaction.mass_weighted = False
    assert interaction.mass_weighted is False


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

    return ParticleInteraction(position, particle_ids, **keywords)


def test_repr(interaction_with_properties):
    expectation = (
        "<ParticleInteraction "
        "position:[0. 0. 0.] "
        "particles:[] "
        "reset_velocities:False "
        "scale:1.0 "
        "mass_weighted:True "
        "max_force:20000.0 "
        "type:gaussian "
        "other:{"
            "'arbitrary_property': 'arbitrary value', "
            "'other_arbitrary_property': 'other arbitrary value'"
        "}>"
    )
    assert repr(interaction_with_properties) == expectation
