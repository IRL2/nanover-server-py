# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Link Narupa's user forces to an OpenMM simulation.

The iMD is hooked to the OpenMM simulation in two places. A
:class:`~simtk.openmm.CustomExternalForce` needs to be added in the system; it
has the components of the iMD force for each atom as a per-atom parameter. A
reporter periodically updates these parameters based on the imd service, and
updates the simulation context.
"""
from typing import Tuple, Dict

import numpy as np

import simtk.openmm as mm
from simtk.openmm import app
from simtk import unit

from narupa.imd.imd_force import calculate_imd_force
from narupa.imd.imd_service import ImdService
from narupa.imd.particle_interaction import ParticleInteraction

NextReport = Tuple[int, bool, bool, bool, bool, bool]


class NarupaImdReporter:
    def __init__(
            self,
            report_interval: int,
            imd_force: mm.CustomExternalForce,
            imd_service: ImdService,
    ):
        self.report_interval = report_interval
        self.imd_force = imd_force
        self.imd_service = imd_service

        # We will not know this values until the beginning of the simulation.
        self.n_particles = None
        self.masses = None
        self.positions = None

        self.is_force_dirty = False

    # The name of the method is part of the OpenMM API. It cannot be made to
    # conform PEP8.
    # noinspection PyPep8Naming
    def describeNextReport(self, simulation: app.Simulation) -> NextReport:
        """
        Called by OpenMM. Indicates when the next report is due and what type
        of data it requires.
        """
        steps = self.report_interval - simulation.currentStep % self.report_interval
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
        if self.masses is None:
            self.n_particles = self.imd_force.getNumParticles()
            self.masses = self.get_masses(simulation.system)
        positions = np.asarray(state.getPositions().value_in_unit(unit.nanometer))
        interactions = self.imd_service.active_interactions
        self._update_forces(positions, interactions, simulation.context)

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


