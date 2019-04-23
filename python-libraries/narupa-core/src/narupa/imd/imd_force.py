"""
Provides a reference implementation of the IMD forces used by Narupa.

"""
from math import exp
from typing import Collection, Tuple

import numpy as np
from narupa.imd.interaction import Interaction


def calculate_imd_force(positions, masses, interactions: Collection[Interaction]) -> Tuple[float, np.array]:
    """
    Reference implementation of the Narupa IMD force.

    :param masses:
    :param positions:
    :param interactions:
    :return: energy in kJ/mol, accumulated forces (in kJ/(mol*nm)) to be applied.
    """

    forces = np.zeros((len(positions), 3))
    total_energy = 0
    for interaction in interactions:
        energy, forces = calculate_single_interaction(positions, masses, interaction, forces)
        total_energy += energy
    return total_energy, forces


def calculate_single_interaction(positions, masses, interaction: Interaction, forces: np.array):
    """
    Calculates the energy and force of a single application of an interaction potential.
    :param masses:
    :param positions:
    :param interaction: An interaction to be applied.
    :param forces: Forces array to accumulate into (in kJ/(mol*nm))
    :return: energy in kJ/mol, accumulated forces (in kJ/(mol*nm)) to be applied.
    """

    centre = get_center_of_mass_subset(positions, masses, interaction.particles)

    try:
        potential_method = interaction_method_map[interaction.type]
    except KeyError:
        if interaction.type is None:
            potential_method = interaction_method_map['gaussian']
        else:
            raise KeyError(f"Unknown interactive force type {interaction.type}.")

    energy, force = potential_method(centre, interaction.position)
    force_per_particle = force / len(interaction.particles)
    energy_per_particle = energy / len(interaction.particles)

    total_energy = 0
    for index in interaction.particles:
        if interaction.mass_weighted:
            mass = masses[index]
        else:
            mass = 1
        total_energy += interaction.scale * mass * energy_per_particle
        forces[index] += interaction.scale * mass * force_per_particle
    return total_energy, forces


def get_center_of_mass_subset(positions, masses, subset=None):
    """
    Gets the center of mass of [a subset of] positions.

    :param positions: List of N vectors representing positions.
    :param masses: List of N vectors representing masses.
    :param subset: Indices [0,N) of positions to include. If None, all positions included.
    :return:
    """
    pos = np.array(positions).reshape((-1, 3))
    if not isinstance(masses, Collection):
        masses = [masses]
    if subset is None:
        subset = range(len(pos))
    com = np.zeros(3)
    total_mass = 0
    for index in subset:
        com += pos[index] * masses[index]
        total_mass += masses[index]
    com = com / total_mass
    return com


def calculate_gaussian_force(particle_position: np.array, interaction_position: np.array, sigma=1):
    """
    \diff{V_{ext}}{\vec{r}_j} = \frac{m_jc}{\sigma^2}\left(\vec{r_j} - \vec{g_i}\right)\exp{\left(\frac{-\left\lVert \vec{r_j} - \vec{g_i} \right\rVert^2}{2\sigma^2}\right)},
    :param particle_position:
    :param interaction_position:
    :param sigma:
    :return:
    """
    # switch to math symbols used in publications.
    r = particle_position
    g = interaction_position
    diff, dist_sqr = _calculate_distances(r, g)
    sigma_sqr = sigma * sigma

    gauss = exp(-dist_sqr / (2 * sigma_sqr))
    energy = -1 * gauss
    # force is negative derivative of energy wrt to position.
    force = diff / sigma_sqr * gauss
    return energy, force


def calculate_spring_force(particle_position: np.array, interaction_position: np.array, k=1):
    r = particle_position
    g = interaction_position

    diff, dist_sqr = _calculate_distances(r, g)
    energy = - k *dist_sqr
    # force is negative derivative of energy wrt to position.
    force = 2 * k * diff
    return energy, force

def _calculate_distances(r, g):
    diff = r - g
    dist_sqr = np.dot(diff, diff)
    return diff, dist_sqr


interaction_method_map = {'gaussian': calculate_gaussian_force, 'spring': calculate_spring_force}
