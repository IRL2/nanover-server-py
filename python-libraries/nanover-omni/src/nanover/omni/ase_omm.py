import warnings
from os import PathLike
from pathlib import Path
from typing import Optional, Callable

import numpy as np
from ase import units, Atoms
from ase.md import Langevin, MDLogger
from ase.md.md import MolecularDynamics
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from openmm.app import Simulation

from nanover.app import NanoverImdApplication
from nanover.ase.converter import (
    EV_TO_KJMOL,
    add_ase_positions_to_frame_data,
    ase_to_frame_data,
)
from nanover.ase.imd_calculator import ImdCalculator
from nanover.ase.openmm import OpenMMCalculator
from nanover.ase.wall_constraint import VelocityWallConstraint
from nanover.omni.ase import InitialState
from nanover.openmm import serializer, openmm_to_frame_data
from nanover.utilities.event import Event


CONSTRAINTS_UNSUPPORTED_MESSAGE = (
    "The simulation contains constraints which will be ignored by this runner!"
)


class ASEOpenMMSimulation:
    """
    A wrapper for ASE OpenMM simulations so they can be run inside the OmniRunner.
    """

    @classmethod
    def from_simulation(
        cls,
        simulation: Simulation,
        *,
        name: Optional[str] = None,
    ):
        """
        Construct this from an existing ASE OpenMM simulation.
        :param simulation: An existing ASE OpenMM Simulation
        :param name: An optional name for the simulation instead of default
        """
        sim = cls(name)
        sim.simulation = simulation
        return sim

    @classmethod
    def from_xml_path(cls, path: PathLike[str], *, name: Optional[str] = None):
        """
        Construct this from an existing NanoVer OpenMM XML file at a given path.
        :param path: Path of the NanoVer OpenMM XML file
        :param name: An optional name for the simulation instead of filename
        """
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
        self.include_velocities = False
        self.include_forces = False
        self.platform: Optional[str] = None

        self.atoms: Optional[Atoms] = None
        self.dynamics: Optional[MolecularDynamics] = None
        self.simulation: Optional[Simulation] = None
        self.openmm_calculator: Optional[OpenMMCalculator] = None
        self.checkpoint: Optional[InitialState] = None

        self.frame_index = 0
        self._frame_adapter: Optional[Callable] = None

    def load(self):
        """
        Load and set up the simulation if it isn't done already.
        """
        if self.simulation is None and self.xml_path is not None:
            with open(self.xml_path) as infile:
                self.simulation = serializer.deserialize_simulation(
                    infile, platform_name=self.platform
                )

        assert self.simulation is not None

        self.openmm_calculator = OpenMMCalculator(self.simulation)
        self.atoms = self.openmm_calculator.generate_atoms()
        if self.use_walls:
            self.atoms.constraints.append(VelocityWallConstraint())
        self.atoms.calc = self.openmm_calculator

        # Set the momenta corresponding to T=300K
        MaxwellBoltzmannDistribution(self.atoms, temperature_K=300)

        self.checkpoint = InitialState(
            positions=self.atoms.get_positions(),
            velocities=self.atoms.get_velocities(),
            cell=self.atoms.get_cell(),
        )

    def reset(self, app_server: NanoverImdApplication):
        """
        Reset the simulation to its initial conditions, reset IMD interactions, and reset frame stream to begin with
        topology and continue.
        :param app_server: The app server hosting the frame publisher and imd state
        """
        assert (
            self.simulation is not None
            and self.atoms is not None
            and self.checkpoint is not None
            and self.openmm_calculator is not None
        )

        self.app_server = app_server

        self.atoms.set_positions(self.checkpoint.positions)
        self.atoms.set_velocities(self.checkpoint.velocities)
        self.atoms.set_cell(self.checkpoint.cell)

        if self.simulation.system.getNumConstraints() > 0:
            warnings.warn(CONSTRAINTS_UNSUPPORTED_MESSAGE)

        # TODO: do we really need to reconstruct the dynamics? (maybe for step count?)
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
            self.openmm_calculator,
            dynamics=self.dynamics,
        )

        # send the initial topology frame
        frame_data = self.make_topology_frame()
        self.app_server.frame_publisher.send_frame(0, frame_data)
        self.frame_index = 1

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

    def advance_by_one_step(self):
        """
        Advance the simulation to the next point a frame should be reported, and send that frame.
        """
        self.advance_to_next_report()

    def advance_by_seconds(self, dt: float):
        """
        Advance playback time by some seconds, and advance the simulation to the next frame output.
        :param dt: Time to advance playback by in seconds (ignored)
        """
        self.advance_to_next_report()

    def advance_to_next_report(self):
        """
        Step the simulation to the next point a frame should be reported, and send that frame.
        """
        assert self.dynamics is not None and self.app_server is not None

        # determine step count for next frame
        steps_to_next_frame = (
            self.frame_interval
            - self.dynamics.get_number_of_steps() % self.frame_interval
        )

        # advance the simulation
        self.dynamics.run(steps_to_next_frame)

        # generate the next frame
        frame_data = self.make_regular_frame()

        # send the next frame
        self.app_server.frame_publisher.send_frame(self.frame_index, frame_data)
        self.frame_index += 1

        # check if excessive energy necessitates reset
        if self.reset_energy is not None and self.app_server is not None:
            energy = self.dynamics.atoms.get_total_energy() * EV_TO_KJMOL
            if not np.isfinite(energy) or energy > self.reset_energy:
                self.on_reset_energy_exceeded.invoke()
                self.reset(self.app_server)

    def make_topology_frame(self):
        """
        Make a NanoVer FrameData corresponding to the current particle positions and topology of the simulation.
        """
        assert self.atoms is not None

        imd_calculator = self.atoms.calc
        topology = imd_calculator.calculator.topology
        frame_data = openmm_to_frame_data(
            state=None,
            topology=topology,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
        )
        add_ase_positions_to_frame_data(frame_data, self.atoms.get_positions())

        return frame_data

    def make_regular_frame(self):
        """
        Make a NanoVer FrameData corresponding to the current state of the simulation.
        """
        assert self.atoms is not None

        frame_data = ase_to_frame_data(self.atoms, topology=False)

        return frame_data
