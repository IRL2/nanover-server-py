# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Link Narupa's user forces to an OpenMM simulation.

The iMD is hooked to the OpenMM simulation in two places. A
:class:`~simtk.openmm.CustomExternalForce` needs to be added in the system; it
has the components of the iMD force for each atom as a per-atom parameter. A
reporter periodically updates these parameters based on the imd service, and
updates the simulation context.

The custom force can be setup using :fun:`create_imd_force` and
:fun:`populate_imd_force`, or using :fun:`add_imd_force_to_system` that combines
the two previous functions. When a simulation is created using
:fun:`narupa.openmm.serializer.deserialize_simulation`, the imd force must be
already present, or must be added by passing it with the ``imd_force``
parameter.

The reporter is :class:`NarupaImdReporter` and both sends the frames and
receives the interactions. It can be use instead of
:class:`narupa.openmm.NarupaReporter` that only sends the frames.

.. code:: python

    from narupa.app import NarupaImdApplication
    from narupa.openmm.serializer import deserialize_simulation
    from narupa.openmm.imd import NarupaImdReporter, create_imd_force

    # Setup the Narupa application server
    # The server is accessible using autoconnect.
    with NarupaImdApplication.basic_server() as app:

        # Create the imd force and a simulation that includes it.
        imd_force = create_imd_force()
        with open('simulation.xml') as infile:
            simulation = deserialize_simulation(infile.read(), imd_force=imd_force)

        # Setup the reporter that does the translation between Narupa and OpenMM
        reporter = NarupaImdReporter(
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
from typing import Tuple, Dict

import numpy as np

import simtk.openmm as mm
from simtk.openmm import app
from simtk import unit

from narupa.imd.imd_force import calculate_imd_force
from narupa.imd.imd_service import ImdService
from narupa.trajectory.frame_publisher import FramePublisher
from narupa.imd.particle_interaction import ParticleInteraction
from .converter import openmm_to_frame_data

NextReport = Tuple[int, bool, bool, bool, bool, bool]


class NarupaImdReporter:
    def __init__(
            self,
            frame_interval: int,
            force_interval: int,
            imd_force: mm.CustomExternalForce,
            imd_service: ImdService,
            frame_publisher: FramePublisher,
    ):
        self.frame_interval = frame_interval
        self.force_interval = force_interval
        self.imd_force = imd_force
        self.imd_service = imd_service
        self.frame_publisher = frame_publisher

        # We will not know these values until the beginning of the simulation.
        self.n_particles = None
        self.masses = None
        self.positions = None

        self.is_force_dirty = False
        self._frame_index = 0

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
        # - not the energies
        # - positions are unwrapped
        return steps, True, False, False, False, True

    def report(self, simulation: app.Simulation, state: mm.State) -> None:
        """
        Called by OpenMM.
        """
        if simulation.currentStep % self.frame_interval == 0:
            frame_data = openmm_to_frame_data(state=state, topology=None)
            self.frame_publisher.send_frame(self._frame_index, frame_data)
            self._frame_index += 1
        if simulation.currentStep % self.force_interval == 0:
            positions = np.asarray(state.getPositions().value_in_unit(unit.nanometer))
            interactions = self.imd_service.active_interactions
            self._update_forces(positions, interactions, simulation.context)

    def _on_first_frame(self, simulation: app.Simulation):
        """
        Do the tasks that are only relevant for the first frame.
        """
        if self.masses is None:
            self.n_particles = self.imd_force.getNumParticles()
            self.masses = self.get_masses(simulation.system)
        if self._frame_index == 0:
            state = simulation.context.getState(getPositions=True)
            topology = simulation.topology
            frame_data = openmm_to_frame_data(state=state, topology=topology)
            self.frame_publisher.send_frame(self._frame_index, frame_data)

    @staticmethod
    def get_masses(system: mm.System) -> np.ndarray:
        """
        Collect the mass, in Dalton, of each particle in an OpenMM system and
        return them as a numpy array.
        """
        return np.array([
            system.getParticleMass(particle).value_in_unit(unit.dalton)
            for particle in range(system.getNumParticles())
        ])

    def _update_forces(
            self,
            positions: np.ndarray,
            interactions: Dict[str, ParticleInteraction],
            context: mm.Context,
    ) -> None:
        """
        Get the forces to apply from the iMD service and communicate them to
        OpenMM.
        """
        context_needs_update = False
        if interactions:
            self._apply_forces(positions, interactions)
            context_needs_update = True
        elif self.is_force_dirty:
            self._reset_forces()
            context_needs_update = True

        if context_needs_update:
            self.imd_force.updateParametersInContext(context)

    def _apply_forces(
            self,
            positions: np.ndarray,
            interactions: Dict[str, ParticleInteraction],
    ):
        """
        Set the iMD forces based on the user interactions.
        """
        _, forces_kjmol = calculate_imd_force(
            positions, self.masses, interactions.values(),
        )
        for particle, force in enumerate(forces_kjmol):
            self.imd_force.setParticleParameters(particle, particle, force)
        self.is_force_dirty = True

    def _reset_forces(self):
        """
        Set all the iMD forces to 0.
        """
        for particle in range(self.n_particles):
            self.imd_force.setParticleParameters(particle, particle, (0, 0, 0))
        self.is_force_dirty = False


def create_imd_force() -> mm.CustomExternalForce:
    """
    Returns an empty OpenMM force to communicate imd forces.

    Each particle in the system has a ``fx``, ``fy``, and ``fz`` parameter to
    provide the arbitrary force components.

    The force needs to be populated to include all the particle in the
    simulation :class:`mm.System`.

    .. seealso: populate_imd_force
    """
    force = mm.CustomExternalForce('-fx * x - fy * y - fz * z')
    force.addPerParticleParameter('fx')
    force.addPerParticleParameter('fy')
    force.addPerParticleParameter('fz')
    return force


def populate_imd_force(force: mm.CustomExternalForce, system: mm.System) -> None:
    """
    Add all the particles to the iMD force.

    The iMD force must be one generated by :func:`create_imd_force`.

    .. seealso: create_imd_force
    """
    for particle in range(system.getNumParticles()):
        force.addParticle(particle, (0, 0, 0))


def add_imd_force_to_system(system: mm.System) -> mm.CustomExternalForce:
    """
    Generate an OpenMM force that accepts arbitrary forces per particle.

    The force is created, populated, added to the system and returned.

    This is the force that is used to communicate the particle interactions from
    Narupa by :class:`NarupaImdReporter`.

    .. seealso: create_imd_force, populate_imd_force
    """
    force = create_imd_force()
    populate_imd_force(force, system)
    system.addForce(force)
    return force
