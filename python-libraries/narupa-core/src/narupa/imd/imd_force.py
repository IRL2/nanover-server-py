# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Provides a reference implementation of the IMD forces used by Narupa.

For details, and if you find these functions helpful, please cite:

.. [1] M. O’Connor et al, “An open-source multi-person virtual reality framework for interactive molecular dynamics:
       from quantum chemistry to drug binding”, arXiv:1902.01827, 2019
"""

from math import exp
from typing import Collection, Tuple

import numpy as np
from narupa.imd.interaction import Interaction


def calculate_imd_force(positions, masses, interactions: Collection[Interaction]) -> Tuple[float, np.array]:
    """
    Reference implementation of the Narupa IMD force.

    Given a collection of interactions, particle positions and masses,
    computes the force to be applied to each particle for each interaction
    and accumulates them into an array.

    :param positions: Collection of N particle positions, in nm.
    :param masses: Collection on N particle masses, in a.m.u
    :param interactions: Collection of interactions to be applied.
    :return: energy in kJ/mol, accumulated forces (in kJ/(mol*nm)) to be applied.
    """

    forces = np.zeros((len(positions), 3))
    total_energy = 0
    for interaction in interactions:
        energy, forces = calculate_single_interaction(positions, masses, interaction, forces)
        total_energy += energy
    return total_energy, forces


def calculate_single_interaction(positions, masses, interaction, forces: np.array) -> Tuple[float, np.array]:
    """
    Calculates the energy and force of a single application of an interaction potential.
    :param positions: Collection of N particle position vectors, in nm.
    :param masses: Collection on N particle masses, in a.m.u
    :param interaction: An interaction to be applied.
    :param forces: Forces array to accumulate into (in kJ/(mol*nm))
    :return: energy in kJ/mol, accumulated forces (in kJ/(mol*nm)) to be applied.
    """

    centre = get_center_of_mass_subset(positions, masses, interaction.particles)

    # fetch the correct potential to use based on the interaction type.
    try:
        potential_method = INTERACTION_METHOD_MAP[interaction.type]
    except KeyError:
        if interaction.type is None:
            potential_method = INTERACTION_METHOD_MAP['gaussian']
        else:
            raise KeyError(f"Unknown interactive force type {interaction.type}.")

    # calculate the overall force to be applied
    energy, force = potential_method(centre, interaction.position)

    # apply to appropriate force to each particle in the selection.
    force_per_particle = force / len(interaction.particles)
    energy_per_particle = energy / len(interaction.particles)
    total_energy, forces = _apply_force_to_particles(forces, energy_per_particle,
                                                     force_per_particle, interaction, masses)
    return total_energy, forces


def _apply_force_to_particles(forces, energy_per_particle, force_per_particle, interaction, masses) \
        -> Tuple[float, np.array]:
    # given the array of forces, energy and force to apply to each particle, applies them, using mass weighting
    # if specified in the interaction.
    total_energy = 0
    for index in interaction.particles:
        if interaction.mass_weighted:
            mass = masses[index]
        else:
            mass = 1
        total_energy += interaction.scale * mass * energy_per_particle
        forces[index] += interaction.scale * mass * force_per_particle
    return total_energy, forces


def get_center_of_mass_subset(positions, masses, subset=None) -> float:
    """
    Gets the center of mass of [a subset of] positions.

    :param positions: List of N vectors representing positions.
    :param masses: List of N vectors representing masses.
    :param subset: Indices [0,N) of positions to include. If None, all positions included.
    :return: The center of mass of the subset of positions.
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


def calculate_gaussian_force(particle_position: np.array, interaction_position: np.array, sigma=1)\
        -> Tuple[float, np.array]:
    """
    Computes the interactive Gaussian force.

    The force applied to the given particle position is determined by the position of a Gaussian centered on the
    interaction position.

    :param particle_position: The position of the particle.
    :param interaction_position: The position of the interaction.
    :param sigma: The width of the Gaussian. Increasing this results in a more diffuse, but longer reaching interaction.
    :return: The energy of the interaction, and the force to be applied to the particle.
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


def calculate_spring_force(particle_position: np.array, interaction_position: np.array, k=1) -> Tuple[float, np.array]:
    """
    Computes the interactive harmonic potential (or spring) force.

    The force applied to the given particle position is determined by placing a spring between the particle position
    and the interaction, and pulling the particle towards the interaction site.

    :param particle_position: The position of the particle.
    :param interaction_position: The position of the interaction.
    :param k: The spring constant. A higher value results in a stronger force.
    :return: The energy of the interaction, and the force to be applied to the particle.
    """
    r = particle_position
    g = interaction_position

    diff, dist_sqr = _calculate_distances(r, g)
    energy = - k * dist_sqr
    # force is negative derivative of energy wrt to position.
    force = 2 * k * diff
    return energy, force

def _calculate_distances(r, g):
    diff = r - g
    dist_sqr = np.dot(diff, diff)
    return diff, dist_sqr

INTERACTION_METHOD_MAP = {'gaussian': calculate_gaussian_force, 'spring': calculate_spring_force}
