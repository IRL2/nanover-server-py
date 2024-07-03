from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Any

from ase import units, Atoms
from ase.md import Langevin
from ase.md.md import MolecularDynamics
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from openmm.app import Simulation

from nanover.app import NanoverImdApplication
from nanover.ase import send_ase_frame
from nanover.ase.imd_calculator import ImdCalculator
from nanover.ase.openmm import OpenMMCalculator
from nanover.ase.wall_constraint import VelocityWallConstraint
from nanover.openmm import serializer


@dataclass
class InitialState:
    positions: Any
    velocities: Any
    cell: Any


class ASEOpenMMSimulation:
    def __init__(self, path: str):
        self.name = Path(path).stem

        self.xml_path = path
        self.app_server: Optional[NanoverImdApplication] = None

        self.frame_interval = 5

        self.atoms: Optional[Atoms] = None
        self.dynamics: Optional[MolecularDynamics] = None
        self.simulation: Optional[Simulation] = None
        self.checkpoint: Optional[InitialState] = None

    def load(self):
        platform = None
        walls = False
        time_step = 1

        with open(str(self.xml_path)) as infile:
            self.simulation = serializer.deserialize_simulation(
                infile.read(), platform_name=platform
            )

        openmm_calculator = OpenMMCalculator(self.simulation)
        self.atoms = openmm_calculator.generate_atoms()
        if walls:
            self.atoms.constraints.append(VelocityWallConstraint())
        self.atoms.calc = openmm_calculator

        # Set the momenta corresponding to T=300K
        MaxwellBoltzmannDistribution(self.atoms, temperature_K=300)

        # We do not remove the center of mass (fixcm=False). If the center of
        # mass translations should be removed, then the removal should be added
        # to the OpenMM system.
        self.dynamics = Langevin(
            atoms=self.atoms,
            timestep=time_step * units.fs,
            temperature_K=300,
            friction=1e-2,
            fixcm=False,
        )

        self.checkpoint = InitialState(
            positions=self.atoms.get_positions(),
            velocities=self.atoms.get_velocities(),
            cell=self.atoms.get_cell(),
        )

    def reset(self, app_server: NanoverImdApplication):
        assert (
            self.dynamics is not None
            and self.atoms is not None
            and self.checkpoint is not None
        )

        self.app_server = app_server
        self.atoms.calc = ImdCalculator(
            self.app_server.imd,
            self.dynamics.atoms.calc,
            dynamics=self.dynamics,
        )

        # replace previous frame method with fresh instance
        self.dynamics.observers.clear()
        self.dynamics.attach(
            send_ase_frame(self.atoms, self.app_server.frame_publisher),
            interval=self.frame_interval,
        )

        self.atoms.set_positions(self.checkpoint.positions)
        self.atoms.set_velocities(self.checkpoint.velocities)
        self.atoms.set_cell(self.checkpoint.cell)

    def advance_by_one_step(self):
        self.advance_to_next_report()

    def advance_by_seconds(self, dt: float):
        self.advance_to_next_report()

    def advance_to_next_report(self):
        assert self.dynamics is not None
        self.dynamics.run(self.frame_interval)
