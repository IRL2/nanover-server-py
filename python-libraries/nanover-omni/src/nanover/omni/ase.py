import warnings
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Optional, Any, Callable

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
from nanover.ase.openmm.runner import openmm_ase_frame_adaptor
from nanover.ase.wall_constraint import VelocityWallConstraint
from nanover.openmm import serializer
from nanover.utilities.event import Event


CONSTRAINTS_UNSUPPORTED_MESSAGE = (
    "The simulation contains constraints which will be ignored by this runner!"
)


@dataclass
class InitialState:
    positions: Any
    velocities: Any
    cell: Any


class ASESimulation:
    """
    A wrapper for ASE simulations so they can be run inside the OmniRunner.
    """

    @classmethod
    def from_dynamics(
        cls,
        dynamics: MolecularDynamics,
        *,
        name: Optional[str] = None,
        frame_method=send_ase_frame,
    ):
        """
        Construct this from an existing ASE dynamics.
        :param dynamics: An existing ASE Dynamics
        :param name: An optional name for the simulation instead of default
        :param frame_method: A
        """
        sim = cls(name)
        sim.dynamics = dynamics
        sim.frame_method = frame_method
        return sim

    @property
    def atoms(self):
        try:
            return self.dynamics.atoms
        except AttributeError:
            return None

    def __init__(self, name: Optional[str] = None):
        self.name = name or "Unnamed ASE OpenMM Simulation"

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

        self.dynamics: Optional[MolecularDynamics] = None
        self.checkpoint: Optional[InitialState] = None

        self._frame_adapter: Optional[Callable] = None

    def load(self):
        """
        Load and set up the simulation if it isn't done already.
        """
        assert self.dynamics is not None

        if self.use_walls:
            self.atoms.constraints.append(VelocityWallConstraint())

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
            self.dynamics is not None
            and self.atoms is not None
            and self.checkpoint is not None
        )

        self.app_server = app_server

        self.atoms.calc = ImdCalculator(
            self.app_server.imd,
            self.atoms.calc,
            dynamics=self.dynamics,
        )

        # remove previous frame adaptor and attach new one
        if self._frame_adapter is not None:
            remove_observer(self.dynamics, self._frame_adapter)

        self._frame_adapter = openmm_ase_frame_adaptor(
            self.atoms,
            self.app_server.frame_publisher,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
        )
        self.dynamics.attach(self._frame_adapter, interval=self.frame_interval)

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

        # reset atoms to initial state
        self.atoms.set_positions(self.checkpoint.positions)
        self.atoms.set_velocities(self.checkpoint.velocities)
        self.atoms.set_cell(self.checkpoint.cell)

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
        Step the simulation to the next point a frame should be reported.
        """
        assert self.dynamics is not None
        self.dynamics.run(self.frame_interval)

        if self.reset_energy is not None and self.app_server is not None:
            energy = self.atoms.get_total_energy() * EV_TO_KJMOL
            if not np.isfinite(energy) or energy > self.reset_energy:
                self.on_reset_energy_exceeded.invoke()
                self.reset(self.app_server)


def remove_observer(dynamics: MolecularDynamics, func: Callable):
    entry = next(entry for entry in dynamics.observers if entry[0] == func)
    try:
        dynamics.observers.remove(entry)
    except StopIteration:
        pass
