"""
Manage an OpenMM CustomExternalForce in conjunction with NanoVer IMD
"""

from collections.abc import Iterable

import numpy as np
import numpy.typing as npt
from openmm import CustomExternalForce, System, unit
from openmm.app import Simulation

from nanover.imd import ImdStateWrapper
from nanover.imd.imd_force import calculate_imd_force, get_sparse_forces
from nanover.imd.particle_interaction import ParticleInteraction
from nanover.trajectory import FrameData

IMD_FORCE_EXPRESSION = "-fx * x - fy * y - fz * z"

ALL_FORCES_GROUP_MASK = 0xFFFFFFFF
IMD_FORCES_GROUP = 31
IMD_FORCES_GROUP_MASK = 1 << IMD_FORCES_GROUP
NON_IMD_FORCES_GROUP_MASK = ALL_FORCES_GROUP_MASK ^ IMD_FORCES_GROUP_MASK


class ImdForceManager:
    def __init__(
        self,
        imd_state: ImdStateWrapper,
        imd_force: CustomExternalForce,
        pbc_vectors: np.ndarray | None,
    ):
        self.imd_state = imd_state
        self.imd_force = imd_force

        self.masses: np.ndarray | None = None
        self.user_forces: np.ndarray = np.empty(self.imd_force.getNumParticles())
        self.total_user_energy = 0.0

        self._prev_particles: set[int] = set()
        self._prev_interactions: dict[str, ParticleInteraction] = {}

        self.periodic_box_lengths: np.ndarray | None = None
        if pbc_vectors is not None:
            # Check that the periodic cell vectors define an orthorhombic cell
            assert np.all(pbc_vectors == np.diagflat(np.diag(pbc_vectors))), (
                "The periodic box vectors do not correspond to an orthorhombic cell. "
                "Periodic boundary conditions are currently only implemented "
                "for orthorhombic systems."
            )
            self.periodic_box_lengths = np.diag(pbc_vectors)

        # clear any residual forces in external force
        for particle in range(self.imd_force.getNumParticles()):
            self.imd_force.setParticleParameters(particle, particle, (0, 0, 0))

    def update_interactions(
        self,
        simulation: Simulation,
        positions: np.ndarray,
    ):
        if self.masses is None:
            self._update_masses(simulation.system)
            assert self.masses is not None

        # which particles and interactions are now active?
        next_interactions = {**self.imd_state.active_interactions}
        next_particles = {
            index
            for interaction in next_interactions.values()
            for index in interaction.particles
        }

        # which previous interactions ended and require velocity reset?
        velocity_resets_interactions = [
            interaction
            for key, interaction in self._prev_interactions.items()
            if key not in next_interactions and interaction.reset_velocities
        ]

        # compute energy and forces for active interactions
        self.total_user_energy, self.user_forces = calculate_imd_force(
            positions,
            self.masses,
            next_interactions.values(),
            self.periodic_box_lengths,
        )

        # overwrite system velocities for velocity reset
        if velocity_resets_interactions:
            state = simulation.context.getState(getVelocities=True)
            velocities = state.getVelocities(asNumpy=True).value_in_unit(
                unit.nanometer / unit.picosecond
            )

            apply_velocity_resets_mean_velocity_removal(
                interactions=velocity_resets_interactions,
                velocities=velocities,
            )

            _ = simulation.integrator.getTemperature()

            # apply_velocity_resets_maxwellboltzmann(
            #     interactions=velocity_resets_interactions,
            #     velocities=velocities,
            #     masses=self.masses,
            #     temperature=simulation.integrator.getTemperature(),
            # )

            simulation.context.setVelocities(velocities)

        # update forces that have changes
        force_resets = self._prev_particles - next_particles
        force_changes = next_particles | force_resets

        if force_changes:
            for particle in force_changes:
                self.imd_force.setParticleParameters(
                    particle, particle, self.user_forces[particle]
                )
            self.imd_force.updateParametersInContext(simulation.context)

        # remember previous interactions and particles
        self._prev_interactions = next_interactions
        self._prev_particles = next_particles

    def add_to_frame_data(self, frame_data: FrameData):
        frame_data.user_energy = self.total_user_energy
        sparse_indices, sparse_forces = get_sparse_forces(self.user_forces)
        frame_data.user_forces_sparse = sparse_forces
        frame_data.user_forces_index = sparse_indices

    def _update_masses(self, system: System):
        self.masses = np.array(
            [
                system.getParticleMass(particle).value_in_unit(unit.dalton)
                for particle in range(system.getNumParticles())
            ]
        )


def apply_velocity_resets_maxwellboltzmann(
    *,
    interactions: Iterable[ParticleInteraction],
    velocities: npt.NDArray,
    masses: npt.NDArray,
    temperature: float,
):
    """Randomise velocities to a certain temperature using a Maxwell-Boltzmann distribution."""
    particles: list[int] = list(
        {particle for interaction in interactions for particle in interaction.particles}
    )

    # TODO: implement
    _ = velocities[particles]
    _ = masses[particles]

    velocities[particles] *= 0


def apply_velocity_resets_mean_velocity_removal(
    *,
    interactions: Iterable[ParticleInteraction],
    velocities: npt.NDArray,
):
    """For each interaction, subtract the mean velocity of those particles."""
    for interaction in interactions:
        mean_velocity = np.average(velocities[interaction.particles], axis=0)
        velocities[interaction.particles] -= mean_velocity


def create_imd_force() -> CustomExternalForce:
    """
    Returns an empty OpenMM force to communicate imd forces.

    Each particle in the system has a ``fx``, ``fy``, and ``fz`` parameter to
    provide the arbitrary force components.

    The force needs to be populated to include all the particle in the
    simulation :class:`mm.System`.

    .. seealso: populate_imd_force
    """
    force = CustomExternalForce(IMD_FORCE_EXPRESSION)
    force.setForceGroup(IMD_FORCES_GROUP)  # Group is used to exclude the force later
    force.addPerParticleParameter("fx")
    force.addPerParticleParameter("fy")
    force.addPerParticleParameter("fz")
    return force


def populate_imd_force(force: CustomExternalForce, system: System) -> None:
    """
    Add all the particles to the iMD force.

    The iMD force must be one generated by :func:`create_imd_force`.

    .. seealso: create_imd_force
    """
    # Attach all the particles to the force object, and set the imd force to 0
    for particle in range(system.getNumParticles()):
        force.addParticle(particle, (0, 0, 0))


def add_imd_force_to_system(system: System) -> CustomExternalForce:
    """
    Generate an OpenMM force that accepts arbitrary forces per particle.

    The force is created, populated, added to the system and returned.

    This is the force that is used to communicate the particle interactions from
    NanoVer by :class:`NanoverImdReporter`.

    .. seealso: create_imd_force, populate_imd_force
    """
    force = create_imd_force()
    populate_imd_force(force, system)
    system.addForce(force)
    return force
