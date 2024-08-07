"""
Link NanoVer's user forces to an OpenMM simulation.

The iMD is hooked to the OpenMM simulation in two places. A
:class:`~openmm.CustomExternalForce` needs to be added in the system; it
has the components of the iMD force for each atom as a per-atom parameter. A
reporter periodically updates these parameters based on the imd service, and
updates the simulation context.

The custom force can be setup using :func:`create_imd_force` and
:func:`populate_imd_force`, or using :func:`add_imd_force_to_system` that combines
the two previous functions. When a simulation is created using
:func:`nanover.openmm.serializer.deserialize_simulation`, the imd force must be
already present, or must be added by passing it with the ``imd_force``
parameter.

The reporter is :class:`NanoverImdReporter` and both sends the frames and
receives the interactions. It can be use instead of
:class:`nanover.openmm.NanoverReporter` that only sends the frames.

.. code:: python

    from nanover.app import NanoverImdApplication
    from nanover.openmm.serializer import deserialize_simulation
    from nanover.openmm.imd import NanoverImdReporter, create_imd_force

    # Setup the NanoVer application server
    # The server is accessible using autoconnect.
    with NanoverImdApplication.basic_server() as app:

        # Create the imd force and a simulation that includes it.
        imd_force = create_imd_force()
        with open('simulation.xml') as infile:
            simulation = deserialize_simulation(infile.read(), imd_force=imd_force)

        # Setup the reporter that does the translation between NanoVer and OpenMM
        reporter = NanoverImdReporter(
            frame_interval=5,
            force_interval=10,
            imd_force=imd_force,
            imd_service=app.imd,
            frame_publisher=app.frame_publisher,
        )
        simulation.reporters.append(reporter)

        # Run the simulation
        while True:
            simulation.run(10)

"""

from typing import Dict, List, Set, Optional, NamedTuple, Tuple
import itertools

import numpy as np
import numpy.typing as npt

from openmm import State, CustomExternalForce, System, Context
from openmm import unit
from openmm.app import Simulation

from nanover.imd.imd_force import calculate_imd_force
from nanover.imd import ImdStateWrapper
from nanover.trajectory.frame_publisher import FramePublisher
from nanover.imd.particle_interaction import ParticleInteraction
from .converter import openmm_to_frame_data
from nanover.trajectory.frame_data import Array2Dfloat, FrameData

IMD_FORCE_EXPRESSION = "-fx * x - fy * y - fz * z"


class NextReport(NamedTuple):
    steps: int
    include_positions: bool
    include_velocities: bool
    include_forces: bool
    include_energies: bool
    wrap_positions: bool


class NanoverImdReporter:
    frame_interval: int
    force_interval: int
    include_velocities: bool
    include_forces: bool
    imd_force: CustomExternalForce
    frame_publisher: FramePublisher
    _frame_index: int

    def __init__(
        self,
        frame_interval: int,
        force_interval: int,
        include_velocities: bool,
        include_forces: bool,
        imd_force: CustomExternalForce,
        imd_state: ImdStateWrapper,
        frame_publisher: FramePublisher,
    ):
        self.frame_interval = frame_interval
        self.force_interval = force_interval
        self.include_velocities = include_velocities
        self.include_forces = include_forces
        self.imd_force = imd_force
        self.frame_publisher = frame_publisher

        self.imd_force_manager = ImdForceManager(imd_state, imd_force)

        self._did_first_frame = False
        self._frame_index = 1

    # The name of the method is part of the OpenMM API. It cannot be made to
    # conform PEP8.
    # noinspection PyPep8Naming
    def describeNextReport(self, simulation: Simulation) -> NextReport:
        """
        Called by OpenMM. Indicates when the next report is due and what type
        of data it requires.
        """
        if not self._did_first_frame:
            self._did_first_frame = True
            self.frame_publisher.send_frame(0, self.make_topology_frame(simulation))

        force_steps = self.force_interval - simulation.currentStep % self.force_interval
        frame_steps = self.frame_interval - simulation.currentStep % self.frame_interval
        steps = min(force_steps, frame_steps)

        return NextReport(
            steps=steps,
            include_positions=True,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
            include_energies=True,
            wrap_positions=False,
        )

    def report(self, simulation: Simulation, state: State) -> None:
        """
        Called by OpenMM.
        """
        positions = state.getPositions(asNumpy=True)
        if simulation.currentStep % self.frame_interval == 0:
            frame_data = self.make_regular_frame(state, positions)
            self.frame_publisher.send_frame(self._frame_index, frame_data)
            self._frame_index += 1
        if simulation.currentStep % self.force_interval == 0:
            self.imd_force_manager.update_interactions(simulation, positions)


    def make_topology_frame(self, simulation: Simulation):
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        topology = simulation.topology
        frame_data = openmm_to_frame_data(state=state, topology=topology)
        return frame_data

    def make_regular_frame(self, state: State, positions: Array2Dfloat):
        frame_data = openmm_to_frame_data(
            state=state,
            topology=None,
            include_positions=False,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
        )
        frame_data.particle_positions = positions
        self.imd_force_manager.add_to_frame_data(frame_data)

        return frame_data


class ImdForceManager:
    def __init__(self, imd_state: ImdStateWrapper, imd_force: CustomExternalForce):
        self.imd_state = imd_state
        self.imd_force = imd_force

        self.masses: Optional[np.ndarray] = None
        self.user_forces: np.ndarray = np.empty(0)
        self.total_user_energy = 0.0

        self._is_force_dirty = False
        self._previous_force_index: Set[int] = set()
        self._total_user_energy = 0.0

    def update_interactions(self, simulation: Simulation, positions: np.ndarray):
        if self.masses is None:
            self._update_masses(simulation.system)

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
    interactions: Dict[str, ParticleInteraction]
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


def get_imd_forces_from_system(system: Simulation) -> List[CustomExternalForce]:
    """
    Find the forces that are compatible with an imd force in a given system.

    A compatible force has the expected energy expression, and contains as
    many particles as the system.

    All the compatible force objects are returned.
    """
    system_num_particles = system.getNumParticles()
    return [
        force
        for force in system.getForces()
        if isinstance(force, CustomExternalForce)
        and force.getEnergyFunction() == IMD_FORCE_EXPRESSION
        and force.getNumParticles() == system_num_particles
    ]


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
