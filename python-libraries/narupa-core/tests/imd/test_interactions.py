from hypothesis import given
from narupa.imd.imd_state import dict_to_interaction, interaction_to_dict
from narupa.imd.particle_interaction import ParticleInteraction

from .. import *


@st.composite
def positions(draw):
    floats = st.floats(allow_nan=False, allow_infinity=False)
    return [draw(floats), draw(floats), draw(floats)]


@st.composite
def particles(draw):
    list = draw(st.lists(st.integers(min_value=0), max_size=5, unique=True))
    list.sort()
    return list


@st.composite
def reset_velocities(draw):
    return draw(st.booleans())


@st.composite
def mass_weighted(draw):
    return draw(st.booleans())


@st.composite
def scale(draw):
    return draw(st.floats(allow_nan=False, allow_infinity=False))


@st.composite
def max_force(draw):
    return draw(st.floats(allow_nan=False))


@st.composite
def interaction_type(draw):
    return draw(st.one_of(st.just('gaussian'), st.just('spring'), st.text()))


@st.composite
def interaction_dictionaries(draw):
    dict = draw(serializable_dictionaries())
    dict['position'] = draw(positions())
    dict['particles'] = draw(particles())
    dict['reset_velocities'] = draw(reset_velocities())
    dict['mass_weighted'] = draw(mass_weighted())
    dict['scale'] = draw(scale())
    dict['max_force'] = draw(max_force())
    dict['interaction_type'] = draw(interaction_type())

    return dict


@st.composite
def interactions(draw):
    dict = draw(serializable_dictionaries())
    return ParticleInteraction(position=draw(positions()),
                               particles=draw(particles()),
                               reset_velocities=draw(reset_velocities()),
                               mass_weighted=draw(mass_weighted()),
                               scale=draw(scale()),
                               max_force=draw(max_force()),
                               interaction_type=draw(interaction_type()),
                               **dict)


@given(interaction_dictionaries())
def test_load_then_save(interaction_dict):
    interaction = dict_to_interaction(interaction_dict)
    new_dict = interaction_to_dict(interaction)
    assert dict_to_interaction(new_dict) == interaction


@given(interactions())
def test_save_then_load(interaction):
    interaction_dict = interaction_to_dict(interaction)
    assert dict_to_interaction(interaction_dict) == interaction


@given(positions(), particles(), reset_velocities(), mass_weighted(), scale(), max_force(), interaction_type(),
       serializable_dictionaries())
def test_constructor(position, particles, reset_velocities, mass_weighted, scale, max_force, interaction_type, extra):
    interaction = ParticleInteraction(position=position,
                                      particles=particles,
                                      reset_velocities=reset_velocities,
                                      mass_weighted=mass_weighted,
                                      scale=scale,
                                      max_force=max_force,
                                      interaction_type=interaction_type,
                                      **extra)
    assert np.allclose(interaction.position, np.array(position))
    assert np.all(interaction.particles == np.array(particles))
    assert interaction.reset_velocities == reset_velocities
    assert interaction.scale == scale
    assert interaction.mass_weighted == mass_weighted
    assert interaction.max_force == max_force
    assert interaction.interaction_type == interaction_type
    for key, value in extra.items():
        assert interaction.properties[key] == value


@given(positions(), positions())
def test_position(orig, new):
    interaction = ParticleInteraction(position=orig)
    assert np.allclose(interaction.position, np.array(orig))
    interaction.position = new
    assert np.allclose(interaction.position, np.array(new))


@given(particles(), particles())
def test_particles(orig, new):
    interaction = ParticleInteraction(particles=orig)
    print(np.array(orig))
    print(interaction.particles)
    assert np.all(interaction.particles == np.array(orig))
    interaction.particles = new
    assert np.all(interaction.particles == np.array(new))


@given(mass_weighted(), mass_weighted())
def test_mass_weighted(orig, new):
    interaction = ParticleInteraction(mass_weighted=orig)
    assert interaction.mass_weighted == orig
    interaction.mass_weighted = new
    assert interaction.mass_weighted == new


@given(reset_velocities(), reset_velocities())
def test_reset_velocities(orig, new):
    interaction = ParticleInteraction(reset_velocities=orig)
    assert interaction.reset_velocities == orig
    interaction.reset_velocities = new
    assert interaction.reset_velocities == new


@given(scale(), scale())
def test_scale(orig, new):
    interaction = ParticleInteraction(scale=orig)
    assert interaction.scale == orig
    interaction.scale = new
    assert interaction.scale == new


@given(max_force(), max_force())
def test_max_force(orig, new):
    interaction = ParticleInteraction(max_force=orig)
    assert interaction.max_force == orig
    interaction.max_force = new
    assert interaction.max_force == new


@given(interaction_type(), interaction_type())
def test_interaction_type(orig, new):
    interaction = ParticleInteraction(interaction_type=orig)
    assert interaction.interaction_type == orig
    interaction.interaction_type = new
    assert interaction.interaction_type == new


@given(serializable_dictionaries(), serializable_dictionaries())
def test_extra(original_extra, new_extra):
    interaction = ParticleInteraction(**original_extra)
    for key, value in original_extra.items():
        assert interaction.properties[key] == value
    for key, value in new_extra.items():
        interaction.properties[key] = value
    for key, value in new_extra.items():
        assert interaction.properties[key] == value


def test_set_particle_unique():
    interaction = ParticleInteraction(particles = [0, 0, 0, 1, 2, 3, 4])
    assert np.all(interaction.particles == [0, 1, 2, 3, 4])
