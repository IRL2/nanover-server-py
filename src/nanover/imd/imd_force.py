"""
Provides a reference implementation of the IMD forces used by NanoVer.

For details, and if you find these functions helpful, please cite [1]_.

.. [1] M. O’Connor et al, “An open-source multi-person virtual reality framework for interactive molecular dynamics:
       from quantum chemistry to drug binding”, arXiv:1902.01827, 2019
"""

import math
from collections.abc import Iterable
from math import exp
from typing import Protocol

import numpy as np
import numpy.typing as npt

from nanover.imd.particle_interaction import ParticleInteraction

_E_SQR = math.sqrt(math.e)


class InvalidInteractionError(ValueError):
    pass


class ForceCalculator(Protocol):
    def __call__(
        self,
        particle_position: npt.NDArray,
        interaction_position: npt.NDArray,
        periodic_box_lengths: npt.NDArray | None = None,
        force_magnitude_limit: float | None = None,
    ) -> tuple[float, npt.NDArray]: ...


def calculate_imd_force(
    positions: npt.NDArray,
    masses: npt.NDArray,
    interactions: Iterable[ParticleInteraction],
    periodic_box_lengths: npt.NDArray | None = None,
) -> tuple[float, npt.NDArray]:
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
    interaction: ParticleInteraction,
    forces: npt.NDArray,
    periodic_box_lengths: npt.NDArray | None = None,
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
        particle_index = interaction.particles[0]
        try:
            center = positions[particle_index]
        except IndexError as e:
            raise InvalidInteractionError(
                f"Particle index {particle_index} out of bounds for {len(positions)} particles."
            ) from e

    # fetch the correct potential to use based on the interaction type.
    try:
        potential_method = INTERACTION_METHOD_MAP[interaction.interaction_type]
    except KeyError:
        raise KeyError(
            f"Unknown interactive force type {interaction.interaction_type}."
        )

    # scale force magnitude limit to account for later scaling
    if np.isfinite(interaction.max_force) and interaction.scale > 0:
        force_magnitude_limit = interaction.max_force / interaction.scale
    else:
        force_magnitude_limit = None

    # calculate the raw (unscaled) force to be applied and associated energy
    raw_energy, raw_force = potential_method(
        particle_position=center,
        interaction_position=interaction.position,
        periodic_box_lengths=periodic_box_lengths,
        force_magnitude_limit=force_magnitude_limit,
    )

    # apply the appropriate force to each particle in the selection.
    total_energy = _apply_force_to_particles(
        forces, raw_energy, raw_force, interaction, masses
    )
    return total_energy


def _apply_force_to_particles(
    forces: np.ndarray,
    raw_energy: float,
    raw_force: np.ndarray,
    interaction: ParticleInteraction,
    masses: np.ndarray,
) -> float:
    """

    Given the array of forces, energy and force to apply to each particle, applies them, using mass weighting
    if specified in the interaction.

    :param forces: array of N particle forces. Interaction force will be added to this array, mutating it.
    :param raw_energy: Raw (unscaled) total interaction energy.
    :param raw_force: Raw (unscaled) total force.
    :param interaction: The interaction being computed.
    :param masses: Array of N masses of the particles.
    :return: The total energy applied.
    """
    particles = interaction.particles
    force_scale = interaction.scale

    if interaction.mass_weighted:
        # distribute weight by particle mass
        interaction_weights = masses[particles]
    else:
        # distribute weight equally over particles with non-zero mass
        interaction_weights = (masses[particles] != 0.0).astype(int)

    total_weight = np.sum(interaction_weights)

    # apply nothing if no particles were weighted
    if total_weight == 0.0:
        interaction_energy = 0.0
        return interaction_energy

    # normalise weights to unit column vector
    interaction_weights = interaction_weights.reshape(-1, 1) / total_weight

    # scale energy by scale factor
    interaction_energy = force_scale * raw_energy
    # scale force and distribute over each particle according to weighting
    interaction_forces = force_scale * raw_force * interaction_weights

    forces[particles] += interaction_forces
    return interaction_energy


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
    periodic_box_lengths: np.ndarray | None = None,
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

    try:
        subset_positions = positions[subset]
    except IndexError as e:
        raise InvalidInteractionError(
            f"Particle indexes {subset} out of bounds for {len(positions)} particles."
        ) from e

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
    periodic_box_lengths: npt.NDArray | None = None,
    force_magnitude_limit: float | None = None,
) -> tuple[float, npt.NDArray]:
    """
    Computes the interactive Gaussian force.

    The force applied to the given particle position is determined by the position of a Gaussian centered on the
    interaction position.

    :param particle_position: The position of the particle.
    :param interaction_position: The position of the interaction.
    :param periodic_box_lengths: Vector of periodic boundary lengths.
    :param force_magnitude_limit: Maximum magnitude permitted for this force.
    :return: The energy of the interaction, and the force to be applied to the particle.
    """
    # The width of the Gaussian. Increasing this results in a more diffuse, but longer reaching interaction.
    sigma = 1

    # vector between particle and interaction, accounting for periodic boundaries
    r = particle_position
    g = interaction_position
    diff, dist_sqr = _calculate_diff_and_sqr_distance(r, g, periodic_box_lengths)

    # energy and force for a gaussian potential
    sigma_sqr = sigma * sigma
    gauss = exp(-dist_sqr / (2 * sigma_sqr))
    energy = 1 - gauss
    # force is negative derivative of energy wrt to position. The minus in the energy cancels with the derivative.
    force = -(diff / sigma_sqr) * gauss

    # scale the entire potential down to limit its peak force
    if force_magnitude_limit is not None:
        # peak force of this potential
        force_magnitude_max = 1 / (sigma * _E_SQR)

        # scale everything down if necessary
        if force_magnitude_max > force_magnitude_limit:
            limit_scale = force_magnitude_limit / force_magnitude_max
            energy *= limit_scale
            force *= limit_scale

    return energy, force


def calculate_spring_force(
    particle_position: npt.NDArray,
    interaction_position: npt.NDArray,
    periodic_box_lengths: npt.NDArray | None = None,
    force_magnitude_limit: float | None = None,
) -> tuple[float, npt.NDArray]:
    """
    Computes the interactive harmonic potential (or spring) force.

    The force applied to the given particle position is determined by placing a spring between the particle position
    and the interaction, and pulling the particle towards the interaction site.

    :param particle_position: The position of the particle.
    :param interaction_position: The position of the interaction.
    :param periodic_box_lengths: Vector of periodic boundary lengths.
    :param force_magnitude_limit: Maximum magnitude permitted for this force.
    :return: The energy of the interaction, and the force to be applied to the particle.
    """
    # The spring constant. A higher value results in a stronger force.
    k = 2

    # vector between particle and interaction, accounting for periodic boundaries
    r = particle_position
    g = interaction_position
    diff, dist_sqr = _calculate_diff_and_sqr_distance(r, g, periodic_box_lengths)

    # limit force by capping the modeled length of the spring
    if force_magnitude_limit is not None:
        # distance at which maximum force is reached
        max_force_distance = force_magnitude_limit / k

        # if distance exceeds max distance, cap distance to max
        if dist_sqr > max_force_distance * max_force_distance:
            diff *= max_force_distance / np.sqrt(dist_sqr)
            dist_sqr = max_force_distance * max_force_distance

    # energy and force for a harmonic potential
    energy = 0.5 * k * dist_sqr
    force = -k * diff

    return energy, force


def calculate_constant_force(
    particle_position: npt.NDArray,
    interaction_position: npt.NDArray,
    periodic_box_lengths: npt.NDArray | None = None,
    force_magnitude_limit: float | None = None,
) -> tuple[float, npt.NDArray]:
    """
    Applies a constant force that is independent of the distance between the particle and the interaction site. Applies
    no force when the two overlap.

    :param particle_position: The position of the particle.
    :param interaction_position: The position of the interaction.
    :param periodic_box_lengths: Vector of periodic boundary lengths.
    :param force_magnitude_limit: Maximum magnitude permitted for this force.
    :return: The energy of the interaction, and the force to be applied to the particle.
    """
    # vector between particle and interaction, accounting for periodic boundaries
    r = particle_position
    g = interaction_position
    diff, dist_sqr = _calculate_diff_and_sqr_distance(r, g, periodic_box_lengths)

    # no energy and force with overlap
    if dist_sqr <= 0:
        force = diff * 0
        energy = 0
        return energy, force

    # force direction
    distance = np.sqrt(dist_sqr)
    unit_force = diff / distance

    # cap force magnitude by limit
    force_magnitude = 1
    if force_magnitude_limit is not None:
        force_magnitude = min(force_magnitude, force_magnitude_limit)

    # energy and force
    energy = float(distance * force_magnitude)
    force = unit_force * force_magnitude

    return energy, force


def _minimum_image(diff, periodic_box_lengths: np.ndarray | None = None) -> np.ndarray:
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
    u: np.ndarray,
    v: np.ndarray,
    periodic_box_lengths: np.ndarray | None = None,
) -> tuple[np.ndarray, float]:
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


def get_sparse_forces(
    user_forces: npt.NDArray,
) -> tuple[npt.NDArray, npt.NDArray]:
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


def expand_sparse_forces(
    atom_count: int,
    sparse_indices: npt.NDArray,
    sparse_forces: npt.NDArray,
) -> npt.NDArray:
    """
    Convert sparse forces of `get_sparse_forces` back into a full array of force per particle.
    """
    user_forces = np.zeros((atom_count, 3), dtype=np.float32)
    for index, force in zip(sparse_indices, sparse_forces):
        user_forces[index, :] = force
    return user_forces


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


INTERACTION_METHOD_MAP: dict[str, ForceCalculator] = {
    "gaussian": calculate_gaussian_force,
    "spring": calculate_spring_force,
    "constant": calculate_constant_force,
}
