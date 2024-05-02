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

from typing import Tuple, Dict, List, Set, Optional
import itertools

import numpy as np
import numpy.typing as npt

import openmm as mm
from openmm import app
from simtk import unit

from nanover.imd.imd_force import calculate_imd_force
from nanover.imd import ImdStateWrapper
from nanover.trajectory.frame_publisher import FramePublisher
from nanover.imd.particle_interaction import ParticleInteraction
from .converter import openmm_to_frame_data

IMD_FORCE_EXPRESSION = "-fx * x - fy * y - fz * z"

NextReport = Tuple[int, bool, bool, bool, bool, bool]


class NanoverImdReporter:
    frame_interval: int
    force_interval: int
    imd_force: mm.CustomExternalForce
    imd_state: ImdStateWrapper
    frame_publisher: FramePublisher
    n_particles: Optional[int]
    masses: Optional[npt.NDArray]
    positions: Optional[npt.NDArray]
    _is_force_dirty: bool
    _previous_force_index: Set[int]
    _frame_index: int

    def __init__(
        self,
        frame_interval: int,
        force_interval: int,
        imd_force: mm.CustomExternalForce,
        imd_state: ImdStateWrapper,
        frame_publisher: FramePublisher,
    ):
        self.frame_interval = frame_interval
        self.force_interval = force_interval
        self.imd_force = imd_force
        self.imd_state = imd_state
        self.frame_publisher = frame_publisher

        # We will not know these values until the beginning of the simulation.
        self.n_particles = None
        self.masses = None
        self.positions = None

        self._is_force_dirty = False
        self._previous_force_index: Set[int] = set()
        self._frame_index = 0
        self._total_user_energy = 0.0

    # The name of the method is part of the OpenMM API. It cannot be made to
    # conform PEP8.
    # noinspection PyPep8Naming
    def describeNextReport(self, simulation: app.Simulation) -> NextReport:
        """
        Called by OpenMM. Indicates when the next report is due and what type
        of data it requires.
        """
        self._on_first_frame(simulation)

        force_steps = self.force_interval - simulation.currentStep % self.force_interval
        frame_steps = self.frame_interval - simulation.currentStep % self.frame_interval
        steps = min(force_steps, frame_steps)
        # The reporter needs:
        # - the positions
        # - not the velocities
        # - not the forces
        # - the energies
        # - positions are unwrapped
        return steps, True, False, False, True, True

    def report(self, simulation: app.Simulation, state: mm.State) -> None:
        """
        Called by OpenMM.
        """
        positions = None
        if simulation.currentStep % self.force_interval == 0:
            interactions = self.imd_state.active_interactions
            positions = state.getPositions(asNumpy=True)
            self._total_user_energy = self._update_forces(
                positions.astype(float), interactions, simulation.context
            )
        if simulation.currentStep % self.frame_interval == 0:
            if positions is None:
                positions = state.getPositions(asNumpy=True)
            try:
                frame_data = openmm_to_frame_data(
                    state=state, topology=None, include_positions=False
                )
            except Exception as err:
                print(f"Error while building a frame: {err}")
            else:
                frame_data.particle_positions = positions
                frame_data.user_energy = self._total_user_energy
                self.frame_publisher.send_frame(self._frame_index, frame_data)
                self._frame_index += 1

    def _on_first_frame(self, simulation: app.Simulation):
        """
        Do the tasks that are only relevant for the first frame.
        """
        if self.masses is None:
            self.n_particles = self.imd_force.getNumParticles()
            self.masses = self.get_masses(simulation.system)
        if self._frame_index == 0:
            state = simulation.context.getState(getPositions=True, getEnergy=True)
            topology = simulation.topology
            try:
                frame_data = openmm_to_frame_data(state=state, topology=topology)
            except Exception as err:
                print(f"Error with the first frame: {err}")
            else:
                self.frame_publisher.send_frame(self._frame_index, frame_data)

    @staticmethod
    def get_masses(system: mm.System) -> np.ndarray:
        """
        Collect the mass, in Dalton, of each particle in an OpenMM system and
        return them as a numpy array.
        """
        return np.array(
            [
                system.getParticleMass(particle).value_in_unit(unit.dalton)
                for particle in range(system.getNumParticles())
            ]
        )

    def _update_forces(
        self,
        positions: np.ndarray,
        interactions: Dict[str, ParticleInteraction],
        context: mm.Context,
    ) -> float:
        """
        Get the forces to apply from the iMD service and communicate them to
        OpenMM.
        """
        energy = 0.0
        context_needs_update = False
        if interactions:
            energy = self._apply_forces(positions, interactions)
            context_needs_update = True
        elif self._is_force_dirty:
            self._reset_forces()
            context_needs_update = True

        if context_needs_update:
            self.imd_force.updateParametersInContext(context)
        return energy

    def _apply_forces(
        self,
        positions: np.ndarray,
        interactions: Dict[str, ParticleInteraction],
    ) -> float:
        """
        Set the iMD forces based on the user interactions.
        """
        if self.masses is None:
            raise InitialisationError
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
        return energy

    def _reset_forces(self):
        """
        Set all the iMD forces to 0.
        """
        for particle in self._previous_force_index:
            self.imd_force.setParticleParameters(particle, particle, (0, 0, 0))
        self._is_force_dirty = False
        self._previous_force_index = set()


class InitialisationError(Exception):
    """
    Error raised when the runner has not been initialised correctly and some
    attribute have not been set.

    This most likely means that `_on_first_frame` has not been called as
    expected.
    """


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


def create_imd_force() -> mm.CustomExternalForce:
    """
    Returns an empty OpenMM force to communicate imd forces.

    Each particle in the system has a ``fx``, ``fy``, and ``fz`` parameter to
    provide the arbitrary force components.

    The force needs to be populated to include all the particle in the
    simulation :class:`mm.System`.

    .. seealso: populate_imd_force
    """
    force = mm.CustomExternalForce(IMD_FORCE_EXPRESSION)
    force.addPerParticleParameter("fx")
    force.addPerParticleParameter("fy")
    force.addPerParticleParameter("fz")
    return force


def populate_imd_force(force: mm.CustomExternalForce, system: mm.System) -> None:
    """
    Add all the particles to the iMD force.

    The iMD force must be one generated by :func:`create_imd_force`.

    .. seealso: create_imd_force
    """
    # Attach all the particles to the force object, and set the imd force to 0
    for particle in range(system.getNumParticles()):
        force.addParticle(particle, (0, 0, 0))


def add_imd_force_to_system(system: mm.System) -> mm.CustomExternalForce:
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


def get_imd_forces_from_system(system: app.Simulation) -> List[mm.CustomExternalForce]:
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
        if isinstance(force, mm.CustomExternalForce)
        and force.getEnergyFunction() == IMD_FORCE_EXPRESSION
        and force.getNumParticles() == system_num_particles
    ]
