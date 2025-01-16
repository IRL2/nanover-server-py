"""
Provides a reference implementation of the IMD forces used by NanoVer.

For details, and if you find these functions helpful, please cite [1]_.

.. [1] M. O’Connor et al, “An open-source multi-person virtual reality framework for interactive molecular dynamics:
       from quantum chemistry to drug binding”, arXiv:1902.01827, 2019
"""

from math import exp
from typing import Tuple, Optional, Iterable, Dict, Protocol

import numpy as np
import numpy.typing as npt
from nanover.imd.particle_interaction import ParticleInteraction


class ForceCalculator(Protocol):
    def __call__(
        self,
        particle_position: npt.NDArray,
        interaction_position: npt.NDArray,
        periodic_box_lengths: Optional[npt.NDArray],
    ) -> Tuple[float, npt.NDArray]: ...


def calculate_imd_force(
    positions: npt.NDArray,
    masses: npt.NDArray,
    interactions: Iterable[ParticleInteraction],
    periodic_box_lengths: Optional[npt.NDArray] = None,
) -> Tuple[float, npt.NDArray]:
    """
    Reference implementation of the NanoVer IMD force.

    Given a collection of interactions, particle positions and masses,
    computes the force to be applied to each particle for each interaction
    and accumulates them into an array.

    :param positions: Array of N particle positions, in nm, with shape (N,3).
    :param masses: Array of N particle masses, in a.m.u, with shape (N,).
    :param interactions: Collection of interactions to be applied.
    :param periodic_box_lengths: Orthorhombic periodic box lengths. If given,
        the minimum image convention is applied to the calculation.
    :return: energy in kJ/mol, accumulated forces (in kJ/(mol*nm)) to be applied.
    """

    forces = np.zeros((len(positions), 3), dtype=np.float32)
    total_energy = 0.0
    for interaction in interactions:
        energy = apply_single_interaction_force(
            positions, masses, interaction, forces, periodic_box_lengths
        )
        total_energy += energy
    return total_energy, forces


def apply_single_interaction_force(
    positions: npt.NDArray,
    masses: npt.NDArray,
    interaction,
    forces: npt.NDArray,
    periodic_box_lengths: Optional[npt.NDArray] = None,
) -> float:
    """
    Calculates the energy and adds the forces to the particles of a single application of an interaction potential.

    :param positions: Collection of N particle position vectors, in nm.
    :param masses: Collection on N particle masses, in a.m.u.
    :param interaction: An interaction to be applied.
    :param forces: Array of N force vectors to accumulate computed forces into (in kJ/(mol*nm)).
    :param periodic_box_lengths: Orthorhombic periodic box lengths to use to apply minimum image convention.
    :return: energy in kJ/mol.
    """

    particle_count = len(interaction.particles)

    if particle_count > 1:
        center = get_center_of_mass_subset(
            positions, masses, interaction.particles, periodic_box_lengths
        )
    else:
        center = positions[interaction.particles[0]]

    # fetch the correct potential to use based on the interaction type.
    try:
        potential_method = INTERACTION_METHOD_MAP[interaction.interaction_type]
    except KeyError:
        raise KeyError(
            f"Unknown interactive force type {interaction.interaction_type}."
        )

    # calculate the overall force to be applied
    energy, force = potential_method(
        center, interaction.position, periodic_box_lengths=periodic_box_lengths
    )
    # apply the appropriate force to each particle in the selection.
    force_per_particle = force / particle_count
    energy_per_particle = energy / particle_count
    total_energy = _apply_force_to_particles(
        forces, energy_per_particle, force_per_particle, interaction, masses
    )
    return total_energy


def _apply_force_to_particles(
    forces: np.ndarray,
    energy_per_particle: float,
    force_per_particle: np.ndarray,
    interaction: ParticleInteraction,
    masses: np.ndarray,
) -> float:
    """

    Given the array of forces, energy and force to apply to each particle, applies them, using mass weighting
    if specified in the interaction.

    :param forces: array of N particle forces. Interaction force will be added to this array, mutating it.
    :param energy_per_particle: Interaction energy per particle.
    :param force_per_particle: Force to apply to each particle.
    :param interaction: The interaction being computed.
    :param masses: Array of N masses of the particles.
    :return: The total energy applied.
    """

    particles = interaction.particles
    scale = interaction.scale
    max_force = interaction.max_force

    if interaction.mass_weighted:
        mass = masses[particles]
        total_mass = mass.sum()
    else:
        mass = np.ones(len(particles))
        total_mass = len(particles)

    total_energy = scale * energy_per_particle * total_mass
    # add the force for each particle, adjusted by mass and scale factor.
    force_to_apply = scale * mass[:, np.newaxis] * force_per_particle[np.newaxis, :]
    # clip the forces into maximum force range.
    force_to_apply_clipped = np.clip(force_to_apply, -max_force, max_force)
    # this is technically incorrect, but deriving the actual energy of a clip will involve a lot of maths
    # for what is essentially just, too much energy.
    total_energy = np.clip(total_energy, -max_force, max_force)
    forces[particles] += force_to_apply_clipped
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


def get_center_of_mass_subset(
    positions: np.ndarray,
    masses: np.ndarray,
    subset: Iterable[int],
    periodic_box_lengths: Optional[np.ndarray] = None,
) -> np.ndarray:
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
    subset = list(subset)
    subset_positions = positions[subset]
    subset_masses = masses[subset, np.newaxis]
    subset_total_mass = subset_masses.sum()
    if subset_total_mass == 0:
        # we raise before actually doing the division since we know it will fail.
        raise ZeroDivisionError(
            "Total mass of subset was zero, cannot compute center of mass!"
        )
    if periodic_box_lengths is not None:
        subset_positions = wrap_pbc(subset_positions, periodic_box_lengths)
    # np.average is slow for small arrays so we use a naive implementation.
    com = (subset_positions * subset_masses).sum(axis=0) / subset_total_mass
    return com


def calculate_gaussian_force(
    particle_position: npt.NDArray,
    interaction_position: npt.NDArray,
    periodic_box_lengths: Optional[npt.NDArray] = None,
) -> Tuple[float, npt.NDArray]:
    """
    Computes the interactive Gaussian force.

    The force applied to the given particle position is determined by the position of a Gaussian centered on the
    interaction position.

    :param particle_position: The position of the particle.
    :param interaction_position: The position of the interaction.
    :param periodic_box_lengths: The periodic box vectors. If passed,
    :return: The energy of the interaction, and the force to be applied to the particle.
    """
    # The width of the Gaussian. Increasing this results in a more diffuse, but longer reaching interaction.
    sigma = 1

    # switch to math symbols used in publications.
    r = particle_position
    g = interaction_position
    diff, dist_sqr = _calculate_diff_and_sqr_distance(r, g, periodic_box_lengths)
    sigma_sqr = sigma * sigma

    gauss = exp(-dist_sqr / (2 * sigma_sqr))
    energy = -gauss
    # force is negative derivative of energy wrt to position. The minus in the energy cancels with the derivative.
    force = -(diff / sigma_sqr) * gauss
    return energy, force


def calculate_spring_force(
    particle_position: npt.NDArray,
    interaction_position: npt.NDArray,
    periodic_box_lengths: Optional[npt.NDArray] = None,
) -> Tuple[float, npt.NDArray]:
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
    # The spring constant. A higher value results in a stronger force.
    k = 2

    r = particle_position
    g = interaction_position

    diff, dist_sqr = _calculate_diff_and_sqr_distance(r, g, periodic_box_lengths)
    energy = 0.5 * k * dist_sqr
    # force is negative derivative of energy wrt to position.
    force = -k * diff
    return energy, force


def calculate_constant_force(
    particle_position: npt.NDArray,
    interaction_position: npt.NDArray,
    periodic_box_lengths: Optional[npt.NDArray] = None,
) -> Tuple[float, npt.NDArray]:
    """
    Applies a constant force that is independent of the distance between the particle and the interaction site. Applies
    no force when the two overlap.

    :param particle_position: The position of the particle.
    :param interaction_position: The position of the interaction.
    :param periodic_box_lengths: Vector of periodic boundary lengths.
    :return: The energy of the interaction, and the force to be applied to the particle.
    """
    distance_vector = _minimum_image(
        interaction_position - particle_position, periodic_box_lengths
    )
    distance_magnitude = np.linalg.norm(distance_vector)

    if distance_magnitude > 0:
        force = distance_vector / distance_magnitude
        energy = 1
    else:
        force = 0
        energy = 0

    return energy, force


def _minimum_image(
    diff, periodic_box_lengths: Optional[np.ndarray] = None
) -> np.ndarray:
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


def _calculate_diff_and_sqr_distance(
    u: np.ndarray, v: np.ndarray, periodic_box_lengths: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, float]:
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
    dist_sqr = float(diff.dot(diff))
    return diff, dist_sqr


def get_sparse_forces(user_forces: npt.NDArray) -> Tuple[npt.NDArray, npt.NDArray]:
    """
    Takes in an array of user forces acting on the system containing N particles
    and outputs two arrays that describe these user forces in a sparse form:

    - The first contains the indices of the particles for which the user forces
      are non-zero
    - The second contains the non-zero user forces associated with each index

    :param user_forces: Array of user forces with dimensions (N, 3)
    :return: Array of particle indices, Array of corresponding user forces
    """
    sparse_indices = np.unique(np.nonzero(user_forces)[0])
    sparse_forces = np.zeros((sparse_indices.shape[0], 3))

    for index in range(sparse_indices.shape[0]):
        sparse_forces[index, :] = user_forces[sparse_indices[index]]

    return sparse_indices, sparse_forces


def calculate_contribution_to_work(forces: npt.NDArray, positions: npt.NDArray):
    r"""
    The expression for the work done on the system by the user is

    .. math::
        W = \sum_{t = 1}^{n_{steps}} \sum_{i = 1}^{N} \mathbf{F}_{i}(t - 1)
         \cdot (\mathbf{r}_{i}(t) - \mathbf{r}_{i}(t - 1)))

    which can be rewritten as

    .. math::
        W = \sum_{t = 1}^{n_{steps}} \bigg(  \sum_{i = 1}^{N} \mathbf{F}_{i}(t - 1)
         \cdot \mathbf{r}_{i}(t) \bigg)  - \bigg(  \sum_{i = 1}^{N} \mathbf{F}_{i}(t - 1)
         \cdot \mathbf{r}_{i}(t - 1) \bigg)

    where the contribution at each value of t is separated into an
    previous-step contribution (t-1) and an on-step contribution (t). Doing so
    enables calculation of the work done on-the-fly without having to save
    the positions of the atoms at each time step that the user applies an
    iMD force.

    This function calculates the contribution to the work done on the system by the user
    for a set of forces and positions, and add it to the work done on the system. Only
    involves the atoms affected by the user interaction.

    :param forces: Array of user forces acting on the system (in NanoVer units of force,
        i.e. kJ mol-1 nm-1)
    :param positions: Array of atomic positions of the atoms on which the user forces
        act (in NanoVer units, i.e. nm)
    :return work_done_contribution: the contribution to the work done on the system (in NanoVer units of energy,
        i.e. kJ mol-1)
    """
    work_done_contribution = 0.0
    for atom in range(len(forces)):
        work_done_contribution += np.dot(np.transpose(forces[atom]), positions[atom])

    return work_done_contribution


INTERACTION_METHOD_MAP: Dict[str, ForceCalculator] = {
    "gaussian": calculate_gaussian_force,
    "spring": calculate_spring_force,
    "constant": calculate_constant_force,
}
