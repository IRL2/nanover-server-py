from pathlib import Path
from queue import Queue
from typing import Optional

from ase import units
from ase.md import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from nanover.app import NanoverImdApplication
from nanover.ase import NanoverASEDynamics, send_ase_frame
from nanover.ase.imd_calculator import ImdCalculator
from nanover.ase.openmm import OpenMMCalculator
from nanover.ase.openmm.runner import openmm_ase_frame_adaptor
from nanover.ase.wall_constraint import VelocityWallConstraint
from nanover.openmm import serializer
from nanover.utilities.timing import VariableIntervalGenerator


class ASEOpenMMSimulation:
    def __init__(self, path: str):
        self.name = Path(path).stem

        self.xml_path = path
        self.app_server: Optional[NanoverImdApplication] = None

        self._variable_interval_generator = VariableIntervalGenerator(1 / 30)

        self.frame_index = 0

        self.frame_interval = 5

    def run(self, app_server: NanoverImdApplication, cancel: Queue):
        platform = None
        simulation_counter = 0
        time_step = 1

        with open(str(self.xml_path)) as infile:
            simulation = serializer.deserialize_simulation(
                infile.read(), platform_name=platform
            )

        walls = False
        openmm_calculator = OpenMMCalculator(simulation)
        atoms = openmm_calculator.generate_atoms()
        if walls:
            atoms.constraints.append(VelocityWallConstraint())
        atoms.calc = openmm_calculator

        # Set the momenta corresponding to T=300K
        MaxwellBoltzmannDistribution(atoms, temperature_K=300)

        # We do not remove the center of mass (fixcm=False). If the center of
        # mass translations should be removed, then the removal should be added
        # to the OpenMM system.
        dynamics = Langevin(
            atoms=atoms,
            timestep=time_step * units.fs,
            temperature_K=300,
            friction=1e-2,
            fixcm=False,
        )

        imd_calculator = ImdCalculator(
            app_server.imd,
            dynamics.atoms.calc,
            dynamics=dynamics,
        )
        atoms.calc = imd_calculator

        # replace previous frame method with fresh instance
        dynamics.attach(
            send_ase_frame(
                atoms, app_server.frame_publisher, simulation_counter
            ),
            interval=self.frame_interval,
        )

        steps: Optional[int] = None
        remaining_steps = steps if steps is not None else float("inf")
        for _ in self._variable_interval_generator.yield_interval():
            if not cancel.empty() or remaining_steps <= 0:
                break
            steps_for_this_iteration = min(self.frame_interval, remaining_steps)
            dynamics.run(steps_for_this_iteration)
            remaining_steps -= steps_for_this_iteration
            #self._reset_if_required(reset_energy)
