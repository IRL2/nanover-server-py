import numpy as np
import pytest
from hypothesis import strategies, given
from math import exp
from narupa.imd.imd_force import (get_center_of_mass_subset, calculate_spring_force,
                                  calculate_gaussian_force, apply_single_interaction_force, calculate_imd_force)
from narupa.imd.particle_interaction import ParticleInteraction

# precomputed results of gaussian force.
EXP_1 = exp(-1 / 2)
EXP_3 = exp(-3 / 2)
UNIT = np.array([1, 1, 1]) / np.linalg.norm([1, 1, 1])


@pytest.fixture
def particle_position():
    return np.array([1, 0, 0])


@pytest.fixture
def interaction_position():
    return np.array([0, 0, 0])


@pytest.fixture
def particles():
    num_particles = 50
    positions = np.array([[i, i, i] for i in range(num_particles)])
    masses = np.array([i + 1 for i in range(num_particles)])
    return positions, masses


@pytest.fixture
def single_interaction():
    position = (0, 0, 0)
    index = 1
    return ParticleInteraction(
        position=position,
        particles=[index],
    )


@pytest.fixture
def single_interaction_multiple_atoms():
    position = (0, 0, 0)
    return ParticleInteraction(
        position=position,
        particles=[1, 2, 3],
    )



@pytest.fixture
def single_interactions():
    num_interactions = 2
    return [single_interaction(position=[i, i, i], index=i) for i in range(num_interactions)]


def test_multiple_interactions(particles):
    """
    Tests multiple concurrent interactions.

    Ensures that equidistant interactions on particles [0,1] and particles [1,2] results in zero force on particle 1,
    and the same (but opposite) forces on atoms 0 and 2.
    """

    positions, masses = particles
    interaction = ParticleInteraction(
        position=[0.5, 0.5, 0.5], particles=[0, 1],
    )
    interaction_2 = ParticleInteraction(
        position=[1.5, 1.5, 1.5], particles=[1, 2],
    )
    # set masses of atoms 0 and 2 to be the same, so things cancel out nicely.
    masses[2] = masses[0]
    single_forces = np.zeros((len(positions), 3))

    single_energy = apply_single_interaction_force(positions, masses, interaction, single_forces)

    energy, forces = calculate_imd_force(positions, masses, [interaction, interaction_2])
    expected_energy = 2 * single_energy
    expected_forces = np.zeros((len(positions), 3))
    expected_forces[0, :] = single_forces[0, :]
    expected_forces[1, :] = 0
    # the uneven masses in this calculation mean both atoms 0 and 2 are pulled towards atom 1.
    expected_forces[2, :] = -single_forces[0, :]

    assert np.allclose(energy, expected_energy)
    assert np.allclose(forces, expected_forces)


@pytest.mark.parametrize("scale", [np.nan, np.infty, -np.infty])
def test_interaction_force_invalid_scale(particles, single_interaction, scale):
    with pytest.raises(ValueError):
        single_interaction.scale = scale


@pytest.mark.parametrize("scale", [-1.0, 0, 100])
def test_interaction_force_single(particles, single_interaction, scale):
    """
    Tests that the interaction force calculation gives the expected result on a single atom, at a particular position,
    with varying scale.
    """
    positions, masses = particles
    forces = np.zeros((len(positions), 3))
    expected_forces = np.zeros((len(positions), 3))
    single_interaction.scale = scale
    energy = apply_single_interaction_force(positions, masses, single_interaction, forces)

    expected_energy = - EXP_3 * scale * masses[single_interaction.particles[0]]
    expected_energy = np.clip(expected_energy,
                              -single_interaction.max_force, single_interaction.max_force)
    expected_forces[1, :] = np.array([-EXP_3 * scale * masses[single_interaction.particles[0]]] * 3)
    expected_forces[1, :] = np.clip(expected_forces[1, :], -single_interaction.max_force,
                                    single_interaction.max_force)

    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(forces, expected_forces, equal_nan=True)


@pytest.mark.parametrize("max_force", [np.nan])
def test_invalid_max_force(single_interaction, max_force):
    with pytest.raises(ValueError):
        single_interaction.max_force = max_force


@pytest.mark.parametrize("max_energy", [0, 1, 1000, np.infty, -np.infty])
def test_interaction_force_max_energy(particles, single_interaction, max_energy):
    """
    Tests that setting the max energy field results in the energy being clamped as expected
    """

    positions, masses = particles
    forces = np.zeros((len(positions), 3))
    expected_forces = np.zeros((len(positions), 3))
    single_interaction.max_force = max_energy
    energy = apply_single_interaction_force(positions, masses, single_interaction, forces)

    expected_energy = - EXP_3 * masses[single_interaction.particles[0]]
    expected_energy = np.clip(expected_energy,
                              -single_interaction.max_force, single_interaction.max_force)
    expected_forces[1, :] = np.array([-EXP_3 * masses[single_interaction.particles[0]]] * 3)
    expected_forces[1, :] = np.clip(expected_forces[1, :], -single_interaction.max_force,
                                    single_interaction.max_force)

    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(forces, expected_forces, equal_nan=True)


@pytest.mark.parametrize("mass", [-1.0, 100, np.nan, np.infty, -np.infty])
def test_interaction_force_mass(particles, single_interaction, mass):
    """
    tests that the interaction force calculation gives the expected result on a single atom, at a particular position,
    with varying mass.
    """
    positions, masses = particles
    forces = np.zeros((len(positions), 3))
    expected_forces = np.zeros((len(positions), 3))
    masses = np.array([mass] * len(masses))
    energy = apply_single_interaction_force(positions, masses, single_interaction, forces)

    expected_energy = np.clip(- EXP_3 * mass, -single_interaction.max_force, single_interaction.max_force)
    expected_forces[1, :] = np.clip(np.array([-EXP_3 * mass] * 3), -single_interaction.max_force,
                                    single_interaction.max_force)

    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(forces, expected_forces, equal_nan=True)


def test_interaction_force_zero_mass_singleatom(particles, single_interaction):
    positions, masses = particles
    forces = np.zeros((len(positions), 3))
    masses = np.array([0.0] * len(masses))

    energy = apply_single_interaction_force(positions, masses, single_interaction, forces)
    assert energy == pytest.approx(0)

def test_interaction_force_zero_mass_multiatom(particles, single_interaction_multiple_atoms):
    positions, masses = particles
    forces = np.zeros((len(positions), 3))
    masses = np.array([0.0] * len(masses))

    with pytest.raises(ZeroDivisionError):
        apply_single_interaction_force(positions, masses, single_interaction_multiple_atoms, forces)


@pytest.mark.parametrize("position, selection, selection_masses",
                         [([1, 1, 1], [0, 1], [1, 2]),
                          ([2, 2, 2], [0, 1], [1, 2]),
                          ([0, 0, 0], [0, 1], [1, 2]),
                          ([0, 0, 0], [0, 1, 49], [1, 2, 10]),
                          ([-5, -5, -5], [0, 1, 49], [1, 2, 10]),
                          ([np.nan, np.nan, np.nan], [0, 1], [1, 2])
                          ])
def test_interaction_force_com(particles, position, selection, selection_masses):
    """
    tests that the interaction force gives the correct result when acting on a group of atoms.
    """
    position = np.array(position)
    selection = np.array(selection)
    interaction = ParticleInteraction(
        position=position,
        particles=selection,)
    positions, masses = particles
    # set non uniform masses based on parameterisation
    for index, mass in zip(selection, selection_masses):
        masses[index] = mass
    forces = np.zeros((len(positions), 3))

    # perform the full calculation to generate expected result.
    com = get_center_of_mass_subset(positions, masses, selection)
    diff = com - interaction.position
    dist_sqr = np.dot(diff, diff)
    expected_energy_per_particle = exp(-dist_sqr / 2) / len(selection)
    expected_energy = sum((-expected_energy_per_particle * masses[index] for index in selection))
    expected_forces = np.zeros((len(positions), 3))
    for index in selection:
        expected_forces[index, :] = -1 * diff * masses[index] * expected_energy_per_particle

    energy = apply_single_interaction_force(positions, masses, interaction, forces)
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(forces, expected_forces, equal_nan=True)


@pytest.mark.parametrize("position, selection, selection_masses",
                         [([1, 1, 1], [0, 1], [1, 2]),
                          ([2, 2, 2], [0, 1], [1, 2]),
                          ([0, 0, 0], [0, 1], [1, 2]),
                          ([0, 0, 0], [0, 1, 49], [1, 2, 10]),
                          ([-5, -5, -5], [0, 1, 49], [1, 2, 10]),
                          ([np.nan, np.nan, np.nan], [0, 1], [1, 2])
                          ])
def test_interaction_force_no_mass_weighting(particles, position, selection, selection_masses):
    """
    tests that the interaction force gives the correct result when acting on a group of atoms.
    """
    position = np.array(position)
    selection = np.array(selection)
    interaction = ParticleInteraction(
        position=position,
        particles=selection,
        mass_weighted=False,
    )
    positions, masses = particles
    # set non uniform masses based on parameterisation
    for index, mass in zip(selection, selection_masses):
        masses[index] = mass
    forces = np.zeros((len(positions), 3))

    # perform the full calculation to generate expected result.
    com = get_center_of_mass_subset(positions, masses, selection)
    diff = com - interaction.position
    dist_sqr = np.dot(diff, diff)
    expected_energy_per_particle = exp(-dist_sqr / 2) / len(selection)
    expected_energy = - sum((expected_energy_per_particle for _ in selection))
    expected_forces = np.zeros((len(positions), 3))
    for index in selection:
        expected_forces[index, :] = -1 * diff * expected_energy_per_particle

    energy = apply_single_interaction_force(positions, masses, interaction, forces)
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(forces, expected_forces, equal_nan=True)


def test_interaction_force_unknown_type(particles, single_interaction):
    single_interaction.interaction_type = "unknown_type"
    positions, masses = particles
    forces = np.zeros((len(positions), 3))

    with pytest.raises(KeyError):
        apply_single_interaction_force(positions, masses, single_interaction, forces)


def test_get_com_all(particles):
    positions, masses = particles
    subset = [i for i in range(len(positions))]

    com = get_center_of_mass_subset(positions, masses, subset)

    expected_com = np.sum(positions * masses[:, None], axis=0) / masses.sum()

    assert np.allclose(com, expected_com)


def test_get_com_subset(particles):
    positions, masses = particles
    subset = [i for i in range(0, len(positions), 2)]

    com = get_center_of_mass_subset(positions, masses, subset)

    expected_com = (
            np.sum(positions[subset] * masses[subset, None], axis=0)
            / masses[subset].sum()
    )

    assert np.allclose(com, expected_com)


@strategies.composite
def random_positions_pbc(draw):
    """
    Generates a random periodic box, and 50 random positions.
    """

    # Generate random polar coordinates and convert them to euclidean
    # coordinates to get a periodic box.
    # box length has to nonzero.
    length = strategies.floats(min_value=0.01, max_value=100,
                               allow_nan=False, allow_infinity=False)

    periodic_box_lengths = np.array([draw(length) for x in range(3)])
    # pick two random points in lowest quadrant of the box.
    # TODO positions at or very near zero cause problems, as the wrap can flip between 0 and box length.
    lengths = np.array([strategies.floats(min_value=0.005, max_value=box_length * 0.5,
                                          allow_nan=False, allow_infinity=False) for box_length in
                        periodic_box_lengths])

    num_particles = 4

    masses = np.array([draw(length) for _ in range(num_particles)])

    positions = np.zeros((num_particles, 3))
    for i in range(num_particles):
        positions[i] = np.array([draw(coord) for coord in lengths])

    # generate random integer values to multiply positions by, putting them in different images.
    images = strategies.integers(min_value=-100, max_value=100)
    image_multiples = np.array([draw(images) for _ in range(3 * num_particles)])
    image_multiples.reshape((num_particles, 3))

    # move points to new random positions around the periodic box.
    positions_periodic = np.zeros((num_particles, 3))
    for i in range(num_particles):
        positions_periodic[i] = positions[i] + image_multiples[i] * periodic_box_lengths

    return positions, masses, positions_periodic, periodic_box_lengths


@given(random_positions_pbc())
def test_get_com_subset_pbc(positions_pbc):
    """
    Tests that the center of mass calculation works when using periodic boundary conditions.
    """
    positions, masses, positions_periodic, periodic_box_lengths = positions_pbc
    subset = [i for i in range(0, len(positions), 2)]

    expected_com = (
            np.sum(positions[subset] * masses[subset, None], axis=0)
            / masses[subset].sum()
    )

    com = get_center_of_mass_subset(positions_periodic, masses, subset, periodic_box_lengths=periodic_box_lengths)

    assert np.allclose(com, expected_com)


def test_get_com_single():
    position = np.array([[1, 0, 0]])
    mass = np.array([20])
    com = get_center_of_mass_subset(position, mass, subset=[0])
    assert np.allclose(com, position)


def test_get_com_different_array_lengths(particles):
    positions, mass = particles
    # make masses array not match positions in length.
    mass = np.array([1])
    subset = range(positions.shape[0])
    with pytest.raises(IndexError):
        get_center_of_mass_subset(positions, mass, subset)


@pytest.mark.parametrize("position, interaction_position, expected_energy, expected_force",
                         [([1, 0, 0], [0, 0, 0], -EXP_1, [-EXP_1, 0, 0]),
                          ([0, 0, 0], [1, 0, 0], -EXP_1, [EXP_1, 0, 0]),
                          ([1, 3, 0], [1, 2, 0], -EXP_1, [0, -EXP_1, 0]),
                          ([1, 3, 3], [1, 3, 2], -EXP_1, [0, 0, -EXP_1]),
                          (UNIT, [0, 0, 0], -EXP_1, np.multiply(UNIT, [-EXP_1, -EXP_1, -EXP_1])),
                          ([1, 2, 3], [1, 2, 3], -1, [0, 0, 0]),
                          ([1, 1, 1], [0, 0, 0], -EXP_3, [-EXP_3] * 3),
                          ([1, 0, 0], [1, 0, 0], -1, [0, 0, 0]),
                          ([-1, -1, -1], [0, 0, 0], -EXP_3, [EXP_3] * 3)])
def test_gaussian_force(position, interaction_position, expected_energy, expected_force):
    # tests gaussian force for various hand evaluated values.
    energy, force = calculate_gaussian_force(np.array(position), np.array(interaction_position))
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(force, expected_force, equal_nan=True)


@pytest.mark.parametrize("position, interaction, expected_energy, expected_force",
                         [([1, 0, 0], [0, 0, 0], 1, [-2, 0, 0]),
                          ([0, 0, 0], [1, 0, 0], 1, [2, 0, 0]),
                          ([1, 3, 0], [1, 2, 0], 1, [0, -2, 0]),
                          ([1, 3, 3], [1, 3, 2], 1, [0, 0, -2]),
                          (UNIT, [0, 0, 0], 1, np.multiply(UNIT, [-2, -2, -2])),
                          ([1, 1, 1], [0, 0, 0], 3, [-2, -2, -2]),
                          ([1, 2, 3], [1, 2, 3], 0, [0, 0, 0]),
                          ([-1, -1, -1], [0, 0, 0], 3, [2, 2, 2])])
def test_spring_force(position, interaction, expected_energy, expected_force):
    energy, force = calculate_spring_force(np.array(position), np.array(interaction))
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(force, expected_force, equal_nan=True)
