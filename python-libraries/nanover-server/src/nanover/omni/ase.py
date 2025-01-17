from dataclasses import dataclass
from typing import Optional, Any, Protocol

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.md import MDLogger
from ase.md.md import MolecularDynamics

from nanover.app import NanoverImdApplication
from nanover.ase.converter import (
    EV_TO_KJMOL,
    ASE_TIME_UNIT_TO_PS,
    ase_atoms_to_frame_data,
    ANG_TO_NM,
    KJMOL_TO_EV,
)
from nanover.ase.imd_calculator import ImdCalculator, ImdForceManager
from nanover.ase.thermo import compute_dof, compute_instantaneous_temperature
from nanover.ase.wall_constraint import VelocityWallConstraint
from nanover.trajectory import FrameData
from nanover.utilities.event import Event
from nanover.imd.imd_force import calculate_contribution_to_work


@dataclass
class InitialState:
    positions: Any
    velocities: Any
    cell: Any


class ASEAtomsToFrameData(Protocol):
    def __call__(self, ase_atoms: Atoms, *, topology: bool, **kwargs) -> FrameData: ...


class ASESimulation:
    """
    A wrapper for ASE simulations so they can be run inside the OmniRunner.
    """

    @classmethod
    def from_ase_dynamics(
        cls,
        dynamics: MolecularDynamics,
        *,
        name: Optional[str] = None,
        ase_atoms_to_frame_data: ASEAtomsToFrameData = ase_atoms_to_frame_data
    ):
        """
        Construct this from an existing ASE dynamics.

        :param dynamics: An existing ASE Dynamics
        :param name: An optional name for the simulation instead of default
        :param ase_atoms_to_frame_data: An optional callback to extra frames from the system
        """
        sim = cls(name)
        sim.dynamics = dynamics
        sim.ase_atoms_to_frame_data = ase_atoms_to_frame_data
        sim.initial_calc = dynamics.atoms.calc

        return sim

    @property
    def atoms(self):
        if self.dynamics is None:
            return None
        else:
            return self.dynamics.atoms

    def __init__(self, name: Optional[str] = None):
        self.name = name or "Unnamed ASE OpenMM Simulation"

        self.app_server: Optional[NanoverImdApplication] = None

        self.on_reset_energy_exceeded = Event()

        self.verbose = False
        self.use_walls = False
        self.reset_energy: Optional[float] = None
        self.frame_interval = 5
        self.include_velocities = False
        self.include_forces = False

        self.dynamics: Optional[MolecularDynamics] = None
        self.checkpoint: Optional[InitialState] = None
        self.initial_calc: Optional[Calculator] = None

        self.imd_calculator: Optional[ImdCalculator] = None

        self.frame_index = 0
        self.ase_atoms_to_frame_data = ase_atoms_to_frame_data

        self.work_done: float = 0.0
        self._work_done_intermediate: float = 0.0
        self._prev_imd_forces: Optional[np.ndarray] = None
        self._prev_imd_indices: Optional[np.ndarray] = None

        self._dof: Optional[int] = None

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
            and self.initial_calc is not None
        )

        self.app_server = app_server

        # reset atoms to initial state
        self.atoms.set_positions(self.checkpoint.positions)
        self.atoms.set_velocities(self.checkpoint.velocities)
        self.atoms.set_cell(self.checkpoint.cell)

        # setup imd calculator
        self.imd_calculator = ImdCalculator(
            self.app_server.imd,
            ImdForceManager(self.app_server.imd, self.atoms),
            self.initial_calc,
            dynamics=self.dynamics,
        )

        # assign imd calculator to atoms
        self.atoms.calc = self.imd_calculator

        self._dof = compute_dof(self.atoms)

        # send the initial topology frame
        frame_data = self.make_topology_frame()
        self.app_server.frame_publisher.send_frame(0, frame_data)
        self.frame_index = 1

        # TODO: deal with this when its clear if dynamics should be reconstructed or not..
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
        assert (
            self.dynamics is not None
            and self.app_server is not None
            and self.imd_calculator is not None
        )

        # determine step count for next frame
        steps_to_next_frame = (
            self.frame_interval
            - self.dynamics.get_number_of_steps() % self.frame_interval
        )

        # advance the simulation
        self.dynamics.run(steps_to_next_frame)

        # Get positions to calculate work done (in NanoVer units)
        positions = self.atoms.get_positions(wrap=False) * ANG_TO_NM

        # Calculate on-step contribution to work
        if self._prev_imd_forces is not None:
            affected_atom_positions = positions[self._prev_imd_indices]
            self._work_done_intermediate += calculate_contribution_to_work(
                self._prev_imd_forces, affected_atom_positions
            )

        # update imd forces and energies
        self.imd_calculator.update_interactions()

        # generate the next frame
        frame_data = self.make_regular_frame()

        # Update work done in frame data
        self.work_done = self._work_done_intermediate
        frame_data.user_work_done = self.work_done

        # Calculate previous-step contribution to work for the next time step
        # (negative contribution, so subtract from the total work done)
        if frame_data.user_forces_sparse is not None:
            affected_atom_positions = positions[frame_data.user_forces_index]
            self._work_done_intermediate -= calculate_contribution_to_work(
                frame_data.user_forces_sparse, affected_atom_positions
            )

        # send the next frame
        self.app_server.frame_publisher.send_frame(self.frame_index, frame_data)
        self.frame_index += 1

        # Update previous step forces (saving them in their sparse form)
        self._prev_imd_forces = frame_data.user_forces_sparse
        self._prev_imd_indices = frame_data.user_forces_index

        # check if excessive energy requires sim reset
        if self.reset_energy is not None and self.app_server is not None:
            energy = self.atoms.get_total_energy() * EV_TO_KJMOL
            if not np.isfinite(energy) or energy > self.reset_energy:
                self.on_reset_energy_exceeded.invoke()
                self.reset(self.app_server)

    def make_topology_frame(self):
        """
        Make a NanoVer FrameData corresponding to the current particle positions and topology of the simulation.
        """
        assert self.atoms is not None

        return self.ase_atoms_to_frame_data(
            self.atoms,
            topology=True,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
        )

    def make_regular_frame(self):
        """
        Make a NanoVer FrameData corresponding to the current state of the simulation.
        """
        assert (
            self.atoms is not None
            and self.imd_calculator is not None
            and self._dof is not None
        )

        frame_data = self.ase_atoms_to_frame_data(
            self.atoms,
            topology=False,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
        )

        # Add simulation time to the frame
        if self.dynamics is not None:
            frame_data.simulation_time = self.dynamics.get_time() * ASE_TIME_UNIT_TO_PS

        # Assume that KE is always available, which is true for this case, and  convert to
        # ASE units for calculating instantaneous temperature
        kinetic_energy = frame_data.kinetic_energy * KJMOL_TO_EV
        frame_data.system_temperature = compute_instantaneous_temperature(
            kinetic_energy, self._dof
        )

        # Add the user forces and user energy to the frame (converting from ASE units
        # to NanoVer units)
        self.imd_calculator.add_to_frame_data(frame_data)

        return frame_data
