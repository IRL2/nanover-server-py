from os import PathLike
from pathlib import Path
from typing import Any

import numpy as np

from openmm.app import Simulation, StateDataReporter

from nanover.app.types import AppServer
from nanover.openmm import serializer, openmm_to_frame_data
from nanover.openmm.imd import (
    create_imd_force,
    add_imd_force_to_system,
    ImdForceManager,
    NON_IMD_FORCES_GROUP_MASK,
)
from nanover.openmm.thermo import compute_instantaneous_temperature, compute_dof
from nanover.trajectory.frame_data import Array2Dfloat
from nanover.imd.imd_force import calculate_contribution_to_work


class OpenMMSimulation:
    """
    A wrapper for OpenMM simulations to run inside the OmniRunner.

    The following attributes can be configured after construction:
    - :attr:`frame_interval`: Number of simulation steps to advance between frames.
    - :attr:`include_velocities`: Include particle velocities in frames.
    - :attr:`include_forces`: Include particle forces in frames.
    - :attr:`platform_name`: Name of OpenMM platform to use when loading the system from XML.
    """

    @classmethod
    def from_simulation(cls, simulation: Simulation, *, name: str | None = None):
        """
        Construct this from an existing OpenMM simulation.

        :param simulation: An existing OpenMM Simulation
        :param name: An optional name for the simulation instead of default
        """
        sim = cls(name)
        sim.simulation = simulation
        sim.imd_force = add_imd_force_to_system(simulation.system)
        sim.simulation.context.reinitialize(preserveState=True)
        sim.get_pbcs()

        sim.checkpoint = sim.simulation.context.createCheckpoint()

        return sim

    @classmethod
    def from_xml_path(cls, path: PathLike[str], *, name: str | None = None):
        """
        Construct this from an existing NanoVer OpenMM XML file at a given path.

        :param path: Path of the NanoVer OpenMM XML file
        :param name: An optional name for the simulation instead of filename
        """
        if name is None:
            name = Path(path).stem
        sim = cls(name)
        sim.xml_path = path
        return sim

    def __init__(self, name: str | None = None):
        self.name = name or "Unnamed OpenMM Simulation"

        self.xml_path: PathLike[str] | None = None
        self.app_server: AppServer | None = None

        self.frame_interval = 5
        """Number of simulation steps to advance between frames."""
        self.include_velocities = False
        """Include particle velocities in frames."""
        self.include_forces = False
        """Include particle forces in frames."""
        self.platform_name: str | None = None
        """Name of OpenMM platform to use at the time the system is loaded from XML."""
        self.use_pbc_wrapping: bool = False
        """Provide atom positions wrapped according to PBC such that each molecule has a center of mass within the
        primary periodic box."""
        self.pbc_vectors: np.ndarray | None = None
        """Array of vectors defining the periodic box used by the simulation (if PBCs are employed)."""

        self.imd_force = create_imd_force()
        self.simulation: Simulation | None = None
        self.checkpoint: Any | None = None
        self.verbose_reporter: StateDataReporter | None = None

        self.frame_index = 0
        self.imd_force_manager: ImdForceManager | None = None

        self.work_done: float = 0.0
        self._work_done_intermediate: float = 0.0
        self._prev_imd_forces: np.ndarray | None = None
        self._prev_imd_indices: np.ndarray | None = None

        self._dof: int | None = None

    def load(self):
        """
        Load and set up the simulation if it isn't done already.
        """
        if self.xml_path is None or self.simulation is not None:
            return

        with open(self.xml_path) as infile:
            self.imd_force = create_imd_force()
            self.simulation = serializer.deserialize_simulation(
                infile, imd_force=self.imd_force, platform_name=self.platform_name
            )

        self.get_pbcs()
        self.checkpoint = self.simulation.context.createCheckpoint()

    def reset(self, app_server: AppServer):
        """
        Reset the simulation to its initial conditions, reset IMD interactions, and reset frame stream to begin with
        topology and continue.

        :param app_server: The app server hosting the frame publisher and imd state
        """
        assert self.simulation is not None and self.checkpoint is not None

        self.app_server = app_server
        self.imd_force_manager = ImdForceManager(self.app_server.imd, self.imd_force)

        self._dof = compute_dof(self.simulation.system)

        # reset imd and work
        self.work_done = 0.0
        self._work_done_intermediate = 0.0
        self._prev_imd_forces = None
        self._prev_imd_indices = None

        # reload initial state and cleanup forces
        self.simulation.context.reinitialize()
        self.simulation.context.loadCheckpoint(self.checkpoint)

        # send the initial topology frame
        frame_data = self.make_topology_frame()
        self.app_server.frame_publisher.send_frame(0, frame_data)
        self.frame_index = 1

        # verbose reporter
        if (
            self.verbose_reporter is not None
            and self.verbose_reporter not in self.simulation.reporters
        ):
            self.simulation.reporters.append(self.verbose_reporter)

    def get_pbcs(self):
        """
        Determine whether the simulation uses periodic boundary conditions and if it does,
        retrieve the periodic box vectors in nanometers.
        """
        assert self.simulation is not None
        self.use_pbc_wrapping = self.simulation.system.usesPeriodicBoundaryConditions()
        if self.use_pbc_wrapping:
            self.pbc_vectors = np.array(
                [
                    [vec3._value.x, vec3._value.y, vec3._value.z]
                    for vec3 in self.simulation.system.getDefaultPeriodicBoxVectors()
                ]
            )

    def advance_by_seconds(self, dt: float):
        """
        Advance playback time by some seconds, and advance the simulation to the next frame output.

        :param dt: Time to advance playback by in seconds (ignored)
        """
        self.advance_to_next_report()

    def advance_by_one_step(self):
        """
        Advance the simulation to the next point a frame should be reported, and send that frame.
        """
        self.advance_to_next_report()

    def advance_to_next_report(self):
        """
        Step the simulation to the next point a frame should be reported, and send that frame.
        """
        assert (
            self.simulation is not None
            and self.imd_force_manager is not None
            and self.app_server is not None
        )

        # determine step count for next frame
        steps_to_next_frame = (
            self.frame_interval - self.simulation.currentStep % self.frame_interval
        )

        # advance the simulation
        self.simulation.step(steps_to_next_frame)

        # fetch positions early, for updating imd
        state = self.simulation.context.getState(
            getPositions=True,
            enforcePeriodicBox=self.use_pbc_wrapping,
        )
        positions = state.getPositions(asNumpy=True)

        # Calculate on-step contribution to work
        if self._prev_imd_forces is not None:
            affected_atom_positions = positions[self._prev_imd_indices]
            self._work_done_intermediate += calculate_contribution_to_work(
                self._prev_imd_forces, affected_atom_positions
            )

        # update imd forces and energies
        self.imd_force_manager.update_interactions(
            self.simulation, positions, self.pbc_vectors
        )

        # generate the next frame with the existing (still valid) positions
        frame_data = self.make_regular_frame(positions)

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

    def make_topology_frame(self):
        """
        Make a NanoVer FrameData corresponding to the current particle positions and topology of the simulation.
        """
        assert self.simulation is not None

        state = self.simulation.context.getState(
            getPositions=True, getEnergy=True, enforcePeriodicBox=self.use_pbc_wrapping
        )
        topology = self.simulation.topology
        frame_data = openmm_to_frame_data(state=state, topology=topology)
        return frame_data

    def make_regular_frame(self, positions: Array2Dfloat | None = None):
        """
        Make a NanoVer FrameData corresponding to the current state of the simulation.

        :param positions: Optionally provided particle positions to save fetching them again.
        """
        assert (
            self.simulation is not None
            and self.imd_force_manager is not None
            and self._dof is not None
        )

        # fetch omm state
        state = self.simulation.context.getState(
            getPositions=positions is None,
            getForces=self.include_forces,
            getVelocities=self.include_velocities,
            getEnergy=True,
            enforcePeriodicBox=self.use_pbc_wrapping,
            groups=NON_IMD_FORCES_GROUP_MASK,
        )

        # generate frame based on basic omm state
        frame_data = openmm_to_frame_data(
            state=state,
            topology=None,
            include_positions=positions is None,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
            state_excludes_imd=True,
        )

        # Assume that the KE is always available, which is true for this case
        frame_data.system_temperature = compute_instantaneous_temperature(
            self.simulation, frame_data.kinetic_energy, self._dof
        )

        # add any provided positions
        if positions is not None:
            frame_data.particle_positions = positions

        # add imd force information
        self.imd_force_manager.add_to_frame_data(frame_data)

        return frame_data
