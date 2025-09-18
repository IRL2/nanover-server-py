"""
Manage an OpenMM CustomExternalForce in conjunction with NanoVer IMD
"""

from typing import Dict, Set, Tuple, Optional
import itertools

import numpy as np
import numpy.typing as npt

from openmm import CustomExternalForce, System, Context
from openmm import unit
from openmm.app import Simulation

from nanover.imd.imd_force import calculate_imd_force, get_sparse_forces
from nanover.imd import ImdStateWrapper
from nanover.imd.particle_interaction import ParticleInteraction
from nanover.trajectory.frame_data import FrameData

IMD_FORCE_EXPRESSION = "-fx * x - fy * y - fz * z"

ALL_FORCES_GROUP_MASK = 0xFFFFFFFF
IMD_FORCES_GROUP = 31
IMD_FORCES_GROUP_MASK = 1 << IMD_FORCES_GROUP
NON_IMD_FORCES_GROUP_MASK = ALL_FORCES_GROUP_MASK ^ IMD_FORCES_GROUP_MASK


class ImdForceManager:
    def __init__(self, imd_state: ImdStateWrapper, imd_force: CustomExternalForce):
        self.imd_state = imd_state
        self.imd_force = imd_force

        self.masses: np.ndarray | None = None
        self.user_forces: np.ndarray = np.empty(0)
        self.total_user_energy = 0.0

        self._is_force_dirty = False
        self._previous_force_index: Set[int] = set()
        self._total_user_energy = 0.0

        self.periodic_box_lengths: np.ndarray | None = None

        # clear any residual forces in external force
        for particle in range(self.imd_force.getNumParticles()):
            self.imd_force.setParticleParameters(particle, particle, (0, 0, 0))

    def update_interactions(
        self,
        simulation: Simulation,
        positions: np.ndarray,
        pbc_vectors: Optional[np.ndarray] = None,
    ):
        if self.masses is None:
            self._update_masses(simulation.system)

        if self.periodic_box_lengths is None and pbc_vectors is not None:
            # Check that the periodic cell vectors define an orthorhombic cell
            assert np.all(pbc_vectors == np.diagflat(np.diag(pbc_vectors))), (
                "The periodic box vectors do not correspond to an orthorhombic cell. "
                "Periodic boundary conditions are currently only implemented "
                "for orthorhombic systems."
            )
            self.periodic_box_lengths = np.diag(pbc_vectors)

        self._update_forces(
            positions.astype(float),
            self.imd_state.active_interactions,
            simulation.context,
        )

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

    def _update_forces(
        self,
        positions: np.ndarray,
        interactions: Dict[str, ParticleInteraction],
        context: Context,
    ) -> Tuple[float, npt.NDArray]:
        """
        Get the forces to apply from the iMD service and communicate them to
        OpenMM.
        """
        energy = 0.0
        forces_kjmol = np.zeros(positions.shape)
        context_needs_update = False
        if interactions:
            energy, forces_kjmol = self._apply_forces(positions, interactions)
            context_needs_update = True
        elif self._is_force_dirty:
            self._reset_forces()
            context_needs_update = True

        if context_needs_update:
            self.imd_force.updateParametersInContext(context)

        self.total_user_energy = energy
        self.user_forces = forces_kjmol

        return energy, forces_kjmol

    def _apply_forces(
        self,
        positions: np.ndarray,
        interactions: Dict[str, ParticleInteraction],
    ) -> Tuple[float, npt.NDArray]:
        """
        Set the iMD forces based on the user interactions.
        """
        assert self.masses is not None
        energy, forces_kjmol = calculate_imd_force(
            positions,
            self.masses,
            interactions.values(),
            self.periodic_box_lengths,
        )
        affected_particles = _build_particle_interaction_index_set(interactions)
        to_reset_particles = self._previous_force_index - affected_particles
        for particle in affected_particles:
            force = forces_kjmol[particle]
            self.imd_force.setParticleParameters(particle, particle, force)
        for particle in to_reset_particles:
            self.imd_force.setParticleParameters(particle, particle, (0, 0, 0))
        self._is_force_dirty = True
        self._previous_force_index = affected_particles
        return energy, forces_kjmol

    def _reset_forces(self):
        """
        Set all the iMD forces to 0.
        """
        for particle in self._previous_force_index:
            self.imd_force.setParticleParameters(particle, particle, (0, 0, 0))
        self._is_force_dirty = False
        self._previous_force_index = set()


def _build_particle_interaction_index_set(
    interactions: Dict[str, ParticleInteraction],
) -> Set[int]:
    """
    Get a set of the indices of the particles involved in interactions.
    """
    indices = (interaction.particles for interaction in interactions.values())
    flatten_indices = itertools.chain(*indices)
    # We need to convert the indices to ints otherwise they are numpy types
    # that protobuf do not support.
    set_of_ints = set(map(int, flatten_indices))
    return set_of_ints


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
