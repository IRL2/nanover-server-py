from math import exp

import numpy as np
import pytest

from narupa.imd.imd_force import _calculate_distances, get_center_of_mass_subset, calculate_spring_force, \
    calculate_gaussian_force, calculate_single_interaction, calculate_imd_force
from narupa.imd.interaction import Interaction

exp_1 = exp(-1 / 2)
exp_3 = exp(-3 / 2)


@pytest.fixture
def particle_position():
    return np.array([1, 0, 0])


@pytest.fixture
def interaction_position():
    return np.array([0, 0, 0])


@pytest.fixture
def particles(num_particles=50):
    positions = np.array([[i, i, i] for i in range(num_particles)])
    masses = np.array([1] * num_particles)
    return positions, masses


@pytest.fixture
def single_interaction(position=(0, 0, 0), index=1):
    return Interaction(position=position, particles=[index])

@pytest.fixture
def single_interactions(num_interactions=2):
    return [single_interaction(position=[i, i, i], index=i) for i in range(num_interactions)]


def test_multiple_interactions(particles):
    """
    Tests multiple concurrent interactions.

    Ensures that equidistant interactions on [0,1] and [1,2] results in zero force on particle 1,
    and expected energy/forces elsewhere.
    :param particles:
    :return:
    """
    interaction = Interaction(position=[0.5,0.5,0.5], particles=[0,1])
    interaction_2 = Interaction(position=[1.5,1.5,1.5], particles=[1,2])
    positions, masses = particles
    single_forces = np.zeros((len(positions),3))
    single_energy, single_forces = calculate_single_interaction(positions, masses, interaction, single_forces)
    energy, forces = calculate_imd_force(positions, masses, [interaction, interaction_2])
    expected_energy = 2 * single_energy
    expected_forces = np.zeros((len(positions), 3))
    expected_forces[0, :] = single_forces[0, :]
    expected_forces[1, :] = 0
    expected_forces[2, :] = single_forces[0, :]
    assert np.allclose(energy, expected_energy)
    assert np.allclose(forces, expected_forces)

@pytest.mark.parametrize("scale", [-1.0, 0, 100, np.nan, np.infty, -np.infty])
def test_interaction_force_single(particles, single_interaction, scale):
    """
    Tests that the interaction force calculation gives the expected result on a single atom, at a particular position,
    with varying scale.
    """
    positions, masses = particles
    forces = np.zeros((len(positions), 3))
    expected_forces = np.zeros((len(positions), 3))
    single_interaction.scale = scale
    energy, forces = calculate_single_interaction(positions, masses, single_interaction, forces)

    expected_energy = -exp_3 * scale
    expected_forces[1, :] = np.array([exp_3 * scale] * 3)

    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(forces, expected_forces, equal_nan=True)

@pytest.mark.parametrize("mass", [-1.0, 0, 100, np.nan, np.infty, -np.infty])
def test_interaction_force_mass(particles, single_interaction, mass):
    """
    tests that the interaction force calculation gives the expected result on a single atom, at a particular position,
    with varying mass.
    """
    positions, masses = particles
    forces = np.zeros((len(positions), 3))
    expected_forces = np.zeros((len(positions), 3))
    masses = np.array([mass] * len(masses))
    energy, forces = calculate_single_interaction(positions, masses, single_interaction, forces)
    # special cases for dividing by zero.
    if mass == 0 or abs(mass) == np.infty:
        expected_energy = np.nan
        expected_forces[1, :] = np.nan
    else:
        expected_energy = -exp_3 * mass
        expected_forces[1, :] = np.array([exp_3 * mass] * 3)
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
def test_interaction_force_com(particles, position, selection, selection_masses):
    """
    tests that the interaction force gives the correct result when acting on a group of atoms.
    """
    position = np.array(position)
    selection = np.array(selection)
    interaction = Interaction(position=position, particles=selection)
    positions, masses = particles
    # set non uniform masses based on parameterisation
    for index, mass in zip(selection, selection_masses):
        masses[index] = mass
    forces = np.zeros((len(positions), 3))

    # perform the full calculation to generate expected result.
    com = get_center_of_mass_subset(positions, masses, selection)
    diff = com - interaction.position
    dist_sqr = np.dot(diff, diff)
    expected_energy_per_particle = - exp(-dist_sqr / 2) / len(selection)
    expected_energy = sum((expected_energy_per_particle * masses[index] for index in selection))
    expected_forces = np.zeros((len(positions), 3))
    for index in selection:
        expected_forces[index, :] = -1 * diff * masses[index] * expected_energy_per_particle

    energy, forces = calculate_single_interaction(positions, masses, interaction, forces)
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
    interaction = Interaction(position=position, particles=selection, mass_weighted=False)
    positions, masses = particles
    # set non uniform masses based on parameterisation
    for index, mass in zip(selection, selection_masses):
        masses[index] = mass
    forces = np.zeros((len(positions), 3))

    # perform the full calculation to generate expected result.
    com = get_center_of_mass_subset(positions, masses, selection)
    diff = com - interaction.position
    dist_sqr = np.dot(diff, diff)
    expected_energy_per_particle = - exp(-dist_sqr / 2) / len(selection)
    expected_energy = sum((expected_energy_per_particle for index in selection))
    expected_forces = np.zeros((len(positions), 3))
    for index in selection:
        expected_forces[index, :] = -1 * diff * expected_energy_per_particle

    energy, forces = calculate_single_interaction(positions, masses, interaction, forces)
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(forces, expected_forces, equal_nan=True)


def test_interaction_force_unknown_type(particles, single_interaction):
    single_interaction.type = "unknown_type"
    positions, masses = particles
    forces = np.zeros((len(positions), 3))

    with pytest.raises(KeyError):
        calculate_single_interaction(positions, masses, single_interaction, forces)


def test_interaction_force_default_type(particles):
    position = [0.5, 0.5, 0.5]
    selection = [0, 1]
    interaction = Interaction(position=position, particles=selection)
    interaction.type = None
    positions, masses = particles
    forces = np.zeros((len(positions), 3))

    expected_forces = np.zeros((len(positions), 3))
    expected_energy = -1
    energy, forces = calculate_single_interaction(positions, masses, interaction, forces)
    assert np.allclose(energy, expected_energy)
    assert np.allclose(forces, expected_forces)

def test_interaction_force_com_zero(particles):
    position = [0.5, 0.5, 0.5]
    selection = [0, 1]
    interaction = Interaction(position=position, particles=selection)
    positions, masses = particles
    forces = np.zeros((len(positions), 3))

    expected_forces = np.zeros((len(positions), 3))
    expected_energy = -1
    energy, forces = calculate_single_interaction(positions, masses, interaction, forces)
    assert np.allclose(energy, expected_energy)
    assert np.allclose(forces, expected_forces)


def test_distance(particle_position, interaction_position):
    diff, dist_sqr = _calculate_distances(particle_position, interaction_position)
    assert dist_sqr == pytest.approx(1.0)
    assert np.allclose(diff, [1, 0, 0])


def test_get_com_all(particles):
    positions, masses = particles
    subset = [i for i in range(len(positions))]
    com = get_center_of_mass_subset(positions, masses, subset)
    assert np.allclose(com, [(len(positions) - 1) / 2] * 3)


def test_get_com_subset(particles):
    positions, masses = particles
    subset = [i for i in range(0, len(positions), 2)]
    subset_positions = [positions[i] for i in subset]
    center = np.mean(subset_positions)
    com = get_center_of_mass_subset(positions, masses, subset)
    assert np.allclose(com, center)


def test_get_com_single():
    position = [1, 0, 0]
    mass = 20
    com = get_center_of_mass_subset(position, mass)
    assert np.allclose(com, position)


def test_get_com_different_lengths(particles):
    positions, mass = particles
    mass = 1
    with pytest.raises(IndexError):
        get_center_of_mass_subset(positions, mass)


@pytest.mark.parametrize("position, interaction, expected_energy, expected_force",
                         [([1, 0, 0], [0, 0, 0], -exp_1, [exp_1, 0, 0]),
                          ([1, 1, 1], [0, 0, 0], -exp_3, [exp_3] * 3),
                          ([1, 0, 0], [1, 0, 0], -1, [0, 0, 0]),
                          ([-1, -1, -1], [0, 0, 0], -exp_3, [-exp_3] * 3)])
def test_gaussian_force(position, interaction, expected_energy, expected_force):
    energy, force = calculate_gaussian_force(np.array(position), np.array(interaction))
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(force, expected_force, equal_nan=True)


@pytest.mark.parametrize("position, interaction, expected_energy, expected_force",
                         [([1, 0, 0], [0, 0, 0], -1, [2, 0, 0]),
                          ([1, 1, 1], [0, 0, 0], -3, [2, 2, 2]),
                          ([1, 0, 0], [1, 0, 0], 0, [0, 0, 0]),
                          ([-1, -1, -1], [0, 0, 0], -3, [-2, -2, -2])])
def test_spring_force(position, interaction, expected_energy, expected_force):
    energy, force = calculate_spring_force(np.array(position), np.array(interaction))
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(force, expected_force, equal_nan=True)
