"""
Provides a reference implementation of the IMD forces used by Narupa.

"""
from math import exp
from typing import Collection, Tuple

import numpy as np

from narupa.ase import converter as converter
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
        energy, forces = calculate_force(positions, masses, interaction, forces)
        total_energy += energy
    return total_energy, forces


def calculate_force(positions, masses, interaction: Interaction, forces: np.array):
    """
    Calculates the energy and force of a single application of an interaction potential.
    :param masses:
    :param positions:
    :param interaction: An interaction to be applied.
    :param forces: Forces array to accumulate into (in kJ/(mol*nm))
    :return: energy in kJ/mol, accumulated forces (in kJ/(mol*nm)) to be applied.
    """

    centre, mass_weighting = get_center_of_mass_subset(positions, masses, interaction.particles)
    if not interaction.mass_weighted:
        mass_weighting = 1

    try:
        potential_method = interaction_method_map[interaction.type]
    except KeyError:
        if interaction.type is None:
            potential_method = interaction_method_map['gaussian']
        else:
            raise KeyError(f"Unknown interactive force type {interaction.type}.")

    energy, force = potential_method(centre, interaction.position, mass=mass_weighting, scale=interaction.scale)
    force_per_particle = force / len(interaction.particles)

    for index in interaction.particles:
        if interaction.mass_weighted:
            mass = masses[index]
        else:
            mass = 1
        forces[index] += mass * force_per_particle
    return energy, forces


def get_center_of_mass_subset(positions, masses, subset):
    """
    Gets the center of mass of a subset of atoms along
    with the total mass.

    :param masses:
    :param positions:
    :param subset:
    :return:
    """

    com = np.zeros(3)
    total_mass = 0
    for index in subset:
        com += positions[index]
        total_mass += masses[index]
    com = com / total_mass
    return com * converter.AngToNm / total_mass, total_mass


def calculate_gaussian_force(particle_position: np.array, interaction_position: np.array, mass=1, scale=1, sigma=1):
    """
    \diff{V_{ext}}{\vec{r}_j} = \frac{m_jc}{\sigma^2}\left(\vec{r_j} - \vec{g_i}\right)\exp{\left(\frac{-\left\lVert \vec{r_j} - \vec{g_i} \right\rVert^2}{2\sigma^2}\right)},
    :param particle_position:
    :param interaction_position:
    :param mass:
    :param scale:
    :param sigma:
    :return:
    """
    # switch to math symbols used in publications.
    r = particle_position
    g = interaction_position
    m = mass
    c = scale
    diff, dist_sqr = _calculate_distances(r, g)
    sigma_sqr = sigma * sigma

    gauss = exp(-dist_sqr / (2 * sigma_sqr))
    energy = -1 * c * m * gauss
    # force is negative derivative of energy wrt to position.
    force = diff / sigma_sqr * c * m * gauss
    return energy, force


def calculate_spring_force(particle_position: np.array, interaction_position: np.array, mass=1, scale=1):
    r = particle_position
    g = interaction_position
    m = mass
    c = scale

    diff, dist_sqr = _calculate_distances(r, g)
    energy = - c * m * dist_sqr
    # force is negative derivative of energy wrt to position.
    force = 2 * c * m * diff
    return energy, force


def _calculate_distances(r, g):
    diff = r - g
    dist_sqr = np.dot(diff, diff)
    return diff, dist_sqr


interaction_method_map = {'gaussian': calculate_gaussian_force, 'spring': calculate_spring_force}
