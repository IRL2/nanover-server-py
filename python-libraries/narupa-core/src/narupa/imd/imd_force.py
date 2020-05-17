# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Provides a reference implementation of the IMD forces used by Narupa.

For details, and if you find these functions helpful, please cite [1]_.

.. [1] M. O’Connor et al, “An open-source multi-person virtual reality framework for interactive molecular dynamics:
       from quantum chemistry to drug binding”, arXiv:1902.01827, 2019
"""

from math import exp
from typing import Iterable, Tuple, Optional

import numpy as np
from narupa.imd.particle_interaction import ParticleInteraction


def calculate_imd_force(positions: np.ndarray, masses: np.ndarray, interactions: Iterable[ParticleInteraction],
                        periodic_box_lengths: Optional[np.ndarray] = None) -> Tuple[float, np.array]:
    """
    Reference implementation of the Narupa IMD force.

    Given a collection of interactions, particle positions and masses,
    computes the force to be applied to each particle for each interaction
    and accumulates them into an array.

    :param positions: Array of N particle positions, in nm, with shape (N,3).
    :param masses: Array of N particle masses, in a.m.u, with shape (N,).
    :param interactions: Collection of interactions to be applied.
    :param periodic_box_lengths: Orthorhombic periodic box lengths. If given, the minimum image convention is applied
        to the calculation.
    :return: energy in kJ/mol, accumulated forces (in kJ/(mol*nm)) to be applied.
    """

    forces = np.zeros((len(positions), 3))
    total_energy = 0
    for interaction in interactions:
        energy = apply_single_interaction_force(positions, masses, interaction, forces, periodic_box_lengths)
        total_energy += energy
    return total_energy, forces


def apply_single_interaction_force(positions: np.ndarray, masses: np.ndarray, interaction, forces: np.ndarray,
                                   periodic_box_lengths: Optional[np.array] = None) -> float:
    """
    Calculates the energy and adds the forces to the particles of a single application of an interaction potential.

    :param positions: Collection of N particle position vectors, in nm.
    :param masses: Collection on N particle masses, in a.m.u
    :param interaction: An interaction to be applied.
    :param forces: Array of N force vectors to accumulate computed forces into (in kJ/(mol*nm))
    :param periodic_box_lengths: Orthorhombic periodic box lengths to use to apply minimum image convention.
    :return: energy in kJ/mol.
    """

    center = get_center_of_mass_subset(positions, masses, interaction.particles, periodic_box_lengths)

    # fetch the correct potential to use based on the interaction type.
    interaction_type = interaction.type if interaction.type is not None else 'gaussian'
    try:
        potential_method = INTERACTION_METHOD_MAP[interaction_type]
    except KeyError:
        raise KeyError(f"Unknown interactive force type {interaction.type}.")

    # calculate the overall force to be applied
    energy, force = potential_method(center, interaction.position, periodic_box_lengths=periodic_box_lengths)
    # apply to appropriate force to each particle in the selection.
    force_per_particle = force / len(interaction.particles)
    energy_per_particle = energy / len(interaction.particles)
    total_energy = _apply_force_to_particles(forces, energy_per_particle,
                                             force_per_particle, interaction, masses)
    return total_energy


def _apply_force_to_particles(forces: np.ndarray,
                              energy_per_particle: float,
                              force_per_particle: np.ndarray,
                              interaction: ParticleInteraction, masses: np.ndarray) -> float:
    """

    Given the array of forces, energy and force to apply to each particle, applies them, using mass weighting
    if specified in the interaction.

    :param forces: array of N particle forces. Interaction force will be added to this array, mutating it.
    :param energy_per_particle: Interaction energy per particle.
    :param force_per_particle: Force to apply to each particle.
    :param interaction: The interaction being computed.
    :param masses: Array of N masses of the particles.
    :return:
    """

    if interaction.mass_weighted:
        mass = masses[interaction.particles]
    else:
        mass = np.ones(len(interaction.particles))
    total_energy = interaction.scale * energy_per_particle * sum(mass)
    # add the force for each particle, adjusted by mass and scale factor.
    force_to_apply = interaction.scale * mass[:, np.newaxis] * force_per_particle[np.newaxis, :]
    # clip the forces into maximum force range.
    force_to_apply_clipped = np.clip(force_to_apply, -interaction.max_force, interaction.max_force)
    # this is technically incorrect, but deriving the actual energy of a clip will involve a lot of maths
    # for what is essentially just, too much energy.
    total_energy = np.clip(total_energy, -interaction.max_force, interaction.max_force)
    forces[interaction.particles] += force_to_apply_clipped
    return total_energy


def wrap_pbc(positions: np.ndarray, periodic_box_lengths: np.ndarray):
    """
    Wraps a list of positions into the given orthorhombic periodic box.

    :param positions: List of N vectors with shape (N,3).
    :param periodic_box_lengths: Box lengths of a periodic box positioned at the origin.
    :return: Positions wrapped into the minimum image of the orthorhombic periodic box.
    """
    # expand the box length vector so it has shape (1,3), so it can be broadcast with the positions.
    box_lengths = periodic_box_lengths[np.newaxis, :]
    wrapped = positions - np.floor(positions / box_lengths) * box_lengths
    return wrapped


def get_center_of_mass_subset(positions: np.ndarray, masses: np.ndarray, subset=None,
                              periodic_box_lengths: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Gets the center of mass of [a subset of] positions.
    If orthorhombic periodic box lengths are given, the minimal image convention is applied, wrapping the subset
    into the periodic boundary before calculating the center of mass.

    :param positions: List of N vectors representing positions.
    :param masses: List of N vectors representing masses.
    :param subset: Indices [0,N) of positions to include. If None, all positions included.
    :param periodic_box_lengths: Orthorhombic periodic box lengths to wrap positions into
        before calculating centre of mass.
    :return: The center of mass of the subset of positions.
    """
    if subset is None:
        subset = range(len(positions))
    pos_subset = positions[subset]
    if periodic_box_lengths is not None:
        pos_subset = wrap_pbc(pos_subset, periodic_box_lengths)
    try:
        com = np.average(pos_subset, weights=masses[subset], axis=0)
    except ZeroDivisionError as e:
        raise ZeroDivisionError("Total mass of subset was zero, cannot compute center of mass!", e)
    return com


def calculate_gaussian_force(particle_position: np.ndarray, interaction_position: np.ndarray, sigma=1,
                             periodic_box_lengths: Optional[np.ndarray] = None) -> Tuple[float, np.ndarray]:
    """
    Computes the interactive Gaussian force.

    The force applied to the given particle position is determined by the position of a Gaussian centered on the
    interaction position.

    :param particle_position: The position of the particle.
    :param interaction_position: The position of the interaction.
    :param periodic_box_lengths: The periodic box vectors. If passed,
    :param sigma: The width of the Gaussian. Increasing this results in a more diffuse, but longer reaching interaction.
    :return: The energy of the interaction, and the force to be applied to the particle.
    """
    # switch to math symbols used in publications.
    r = particle_position
    g = interaction_position
    diff, dist_sqr = _calculate_diff_and_sqr_distance(r, g, periodic_box_lengths)
    sigma_sqr = sigma * sigma

    gauss = exp(-dist_sqr / (2 * sigma_sqr))
    energy = - gauss
    # force is negative derivative of energy wrt to position. the minus in the energy cancels with the derivative.
    force = - (diff / sigma_sqr) * gauss
    return energy, force


def calculate_spring_force(particle_position: np.array, interaction_position: np.array, k=1,
                           periodic_box_lengths: Optional[np.ndarray] = None) -> Tuple[float, np.array]:
    """
    Computes the interactive harmonic potential (or spring) force.

    The force applied to the given particle position is determined by placing a spring between the particle position
    and the interaction, and pulling the particle towards the interaction site.

    :param particle_position: The position of the particle.
    :param interaction_position: The position of the interaction.
    :param k: The spring constant. A higher value results in a stronger force.
    :param periodic_box_lengths: Vector of periodic boundary lengths.
    :return: The energy of the interaction, and the force to be applied to the particle.
    """
    r = particle_position
    g = interaction_position

    diff, dist_sqr = _calculate_diff_and_sqr_distance(r, g, periodic_box_lengths)
    energy = k * dist_sqr
    # force is negative derivative of energy wrt to position.
    force = - 2 * k * diff
    return energy, force


def _minimum_image(diff, periodic_box_lengths: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Gets the difference between two vectors under minimum image convention for a cubic periodic box.
    :diff The difference between two vectors.
    :param periodic_box_lengths: Vector of length 3 of box lengths for an orthorhombic periodic boundary.
    :return:
    """
    if periodic_box_lengths is not None:
        pbc_recipricol = np.reciprocal(periodic_box_lengths)
        rounded = np.round(diff * pbc_recipricol)
        diff -= periodic_box_lengths * rounded
    return diff


def _calculate_diff_and_sqr_distance(u: np.ndarray, v: np.ndarray,
                                     periodic_box_lengths: Optional[np.ndarray] = None) -> Tuple[np.ndarray, float]:
    """
    Calculates the difference and square of the distance between two vectors.
    A utility function for computing gradients based on this distance.
    :param u: Vector of length N.
    :param v: Vector of length N.
    :param periodic_box_lengths: Vector of length 3 of box lengths for an orthorhombic periodic boundary. If passed,
    minimum image convention will be used.
    :return: Tuple consisting of the difference between r and g and the square magnitude between them.
    """
    diff = u - v
    diff = _minimum_image(diff, periodic_box_lengths)
    dist_sqr = float(np.dot(diff, diff))
    return diff, dist_sqr


INTERACTION_METHOD_MAP = {'gaussian': calculate_gaussian_force, 'spring': calculate_spring_force}
