from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Optional, Any

import numpy as np
from ase import units, Atoms
from ase.md import Langevin, MDLogger
from ase.md.md import MolecularDynamics
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from openmm.app import Simulation

from nanover.app import NanoverImdApplication
from nanover.ase import send_ase_frame
from nanover.ase.converter import EV_TO_KJMOL
from nanover.ase.imd_calculator import ImdCalculator
from nanover.ase.openmm import OpenMMCalculator
from nanover.ase.wall_constraint import VelocityWallConstraint
from nanover.openmm import serializer
from nanover.utilities.event import Event


@dataclass
class InitialState:
    positions: Any
    velocities: Any
    cell: Any


class ASEOpenMMSimulation:
    @classmethod
    def from_simulation(cls, simulation: Simulation, *, name: Optional[str] = None):
        sim = cls(name)
        sim.simulation = simulation
        return sim

    @classmethod
    def from_xml_path(cls, path: PathLike[str], *, name: Optional[str] = None):
        sim = cls(name or Path(path).stem)
        sim.xml_path = path
        return sim

    def __init__(self, name: Optional[str] = None):
        self.name = name or "Unnamed ASE OpenMM Simulation"

        self.xml_path: Optional[PathLike[str]] = None
        self.app_server: Optional[NanoverImdApplication] = None

        self.on_reset_energy_exceeded = Event()

        self.verbose = False
        self.use_walls = False
        self.reset_energy: Optional[float] = None
        self.time_step = 1
        self.frame_interval = 5
        self.platform: Optional[str] = None

        self.atoms: Optional[Atoms] = None
        self.dynamics: Optional[MolecularDynamics] = None
        self.simulation: Optional[Simulation] = None
        self.checkpoint: Optional[InitialState] = None

    def load(self):
        if self.simulation is None and self.xml_path is not None:
            with open(self.xml_path) as infile:
                self.simulation = serializer.deserialize_simulation(
                    infile, platform_name=self.platform
                )

        assert self.simulation is not None

        openmm_calculator = OpenMMCalculator(self.simulation)
        self.atoms = openmm_calculator.generate_atoms()
        if self.use_walls:
            self.atoms.constraints.append(VelocityWallConstraint())
        self.atoms.calc = openmm_calculator

        # Set the momenta corresponding to T=300K
        MaxwellBoltzmannDistribution(self.atoms, temperature_K=300)

        self.checkpoint = InitialState(
            positions=self.atoms.get_positions(),
            velocities=self.atoms.get_velocities(),
            cell=self.atoms.get_cell(),
        )

    def reset(self, app_server: NanoverImdApplication):
        assert (
            self.atoms is not None
            and self.checkpoint is not None
        )

        self.app_server = app_server

        # We do not remove the center of mass (fixcm=False). If the center of
        # mass translations should be removed, then the removal should be added
        # to the OpenMM system.
        self.dynamics = Langevin(
            atoms=self.atoms,
            timestep=self.time_step * units.fs,
            temperature_K=300,
            friction=1e-2,
            fixcm=False,
        )

        self.atoms.calc = ImdCalculator(
            self.app_server.imd,
            self.dynamics.atoms.calc,
            dynamics=self.dynamics,
        )

        self.dynamics.attach(
            send_ase_frame(self.atoms, self.app_server.frame_publisher),
            interval=self.frame_interval,
        )

        if self.verbose:
            self.dynamics.attach(
                MDLogger(
                    self.dynamics,
                    self.atoms,
                    "-",
                    header=True,
                    stress=False,
                    peratom=False,
                ),
                interval=100,
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

        if self.reset_energy is not None and self.app_server is not None:
            energy = self.dynamics.atoms.get_total_energy() * EV_TO_KJMOL
            if not np.isfinite(energy) or energy > self.reset_energy:
                self.on_reset_energy_exceeded.invoke()
                self.reset(self.app_server)
