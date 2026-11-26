from collections.abc import Callable
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any

import numpy as np
from openmm import CustomExternalForce
from openmm.app import Simulation
from openmm.unit import nanometer

from nanover.core import AppServer, Simulation as NanoverSimulation
from nanover.imd.imd_force import calculate_contribution_to_work
from nanover.utilities.event import Event

from . import serializer
from .bundles import unbundle_openmm_simulation
from .converter import openmm_to_frame_data
from .imd import (
    NON_IMD_FORCES_GROUP_MASK,
    ImdForceManager,
    add_imd_force_to_system,
)
from .thermo import compute_dof, compute_instantaneous_temperature


class OpenMMSimulation(NanoverSimulation):
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
        loader = Loader[OpenMMSimulation](
            lambda _: LoadedSetup.from_simulation(simulation)
        )
        return cls(loader, name)

    @classmethod
    def from_bundle_path(cls, path: str | PathLike[str], *, name: str | None = None):
        loader = Loader[OpenMMSimulation](lambda _: LoadedSetup.from_bundle_path(path))
        return cls(loader, name if name is not None else Path(path).stem)

    @classmethod
    def from_xml_path(cls, path: str | PathLike[str], *, name: str | None = None):
        """
        Construct this from an existing NanoVer OpenMM XML file at a given path.

        :param path: Path of the NanoVer OpenMM XML file
        :param name: An optional name for the simulation instead of filename
        """
        loader = Loader[OpenMMSimulation](lambda _: LoadedSetup.from_bundle_path(path))
        return cls(loader, name if name is not None else Path(path).stem)

    def __init__(self, loader: "Loader[OpenMMSimulation]", name: str | None = None):
        self.name = name or "Unnamed OpenMM Simulation"

        self.xml_path: PathLike[str] | None = None
        self.app_server: AppServer | None = None

        self.loader = loader

        self.frame_interval = 5
        """Number of simulation steps to advance between frames."""
        self.include_velocities = False
        """Include particle velocities in frames."""
        self.include_forces = False
        """Include particle forces in frames."""
        self.platform_name: str | None = None
        """Name of OpenMM platform to use at the time the system is loaded from XML."""
        self.use_pbc_wrapping: bool | None = None
        """Provide atom positions wrapped according to PBC such that each molecule has a center of mass within the
        primary periodic box."""
        self.reset_energy: float | None = None
        """System's kinetic energy threshold for automatic reset"""
        self.on_reset_energy_exceeded = Event()
        self.imd_force_manager: ImdForceManager | None = None

        self.work_done: float = 0.0
        self._work_done_intermediate: float = 0.0
        self._prev_imd_forces: np.ndarray | None = None
        self._prev_imd_indices: np.ndarray | None = None

    def load(self):
        """
        Load and set up the simulation if it isn't done already.
        """
        return self.loader.load(self)

    def reset(self, app_server: AppServer):
        """
        Reset the simulation to its initial conditions, reset IMD interactions, and reset frame stream to begin with
        topology and continue.

        :param app_server: The app server hosting the frame publisher and imd state
        """
        setup = self.loader.load(self)

        self.app_server = app_server
        self.imd_force_manager = ImdForceManager(
            self.app_server.imd,
            setup.imd_force,
            setup.pbc_vectors if self.use_pbc_wrapping else None,
        )

        # reset imd and work
        self.work_done = 0.0
        self._work_done_intermediate = 0.0
        self._prev_imd_forces = None
        self._prev_imd_indices = None

        # reload initial state and cleanup forces
        setup.simulation.context.reinitialize()
        setup.simulation.context.loadCheckpoint(setup.checkpoint)

        # send the initial topology frame
        frame_data = self.make_topology_frame()
        frame_data.simulation_name = self.name
        app_server.frame_publisher.send_clear()
        app_server.frame_publisher.send_frame(frame_data)

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
        assert self.imd_force_manager is not None and self.app_server is not None

        setup = self.load()

        # determine step count for next frame
        steps_to_next_frame = (
            self.frame_interval - setup.simulation.currentStep % self.frame_interval
        )

        # advance the simulation
        setup.simulation.step(steps_to_next_frame)

        # fetch positions early, for updating imd
        state = setup.simulation.context.getState(
            getPositions=True,
            enforcePeriodicBox=self.use_pbc_wrapping or False,
        )
        positions = state.getPositions(asNumpy=True).value_in_unit(nanometer)

        # Calculate on-step contribution to work
        if self._prev_imd_forces is not None:
            affected_atom_positions = positions[self._prev_imd_indices]
            self._work_done_intermediate += calculate_contribution_to_work(
                self._prev_imd_forces, affected_atom_positions
            )

        # update imd forces and energies
        self.imd_force_manager.update_interactions(setup.simulation, positions)

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
        self.app_server.frame_publisher.send_frame(frame_data)

        # Update previous step forces (saving them in their sparse form)
        self._prev_imd_forces = frame_data.user_forces_sparse
        self._prev_imd_indices = frame_data.user_forces_index

        # check if excessive energy requires sim reset
        if self.reset_energy is not None and self.app_server is not None:
            energy = frame_data.kinetic_energy
            if not np.isfinite(energy) or energy > self.reset_energy:
                self.on_reset_energy_exceeded.invoke()
                self.reset(self.app_server)

    def make_topology_frame(self):
        """
        Make a NanoVer FrameData corresponding to the current particle positions and topology of the simulation.
        """
        setup = self.load()

        state = setup.simulation.context.getState(
            getPositions=True,
            getEnergy=True,
            enforcePeriodicBox=self.use_pbc_wrapping or False,
        )
        topology = setup.simulation.topology
        frame_data = openmm_to_frame_data(state=state, topology=topology)
        return frame_data

    def make_regular_frame(self, positions: np.ndarray | None = None):
        """
        Make a NanoVer FrameData corresponding to the current state of the simulation.

        :param positions: Optionally provided particle positions to save fetching them again.
        """
        assert self.imd_force_manager is not None

        setup = self.load()

        # fetch omm state
        state = setup.simulation.context.getState(
            getPositions=positions is None,
            getForces=self.include_forces,
            getVelocities=self.include_velocities,
            getEnergy=True,
            enforcePeriodicBox=self.use_pbc_wrapping or False,
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
            setup.simulation, frame_data.kinetic_energy, setup.dof
        )

        # add any provided positions
        if positions is not None:
            frame_data.particle_positions = positions.astype(np.float32)

        # add imd force information
        self.imd_force_manager.add_to_frame_data(frame_data)

        return frame_data


@dataclass(kw_only=True)
class LoadedSetup:
    simulation: Simulation
    checkpoint: Any
    imd_force: CustomExternalForce
    pbc_vectors: np.ndarray | None
    dof: int

    @classmethod
    def from_bundle_path(cls, path: str | PathLike[str]):
        return LoadedSetup.from_simulation(unbundle_openmm_simulation(path))

    @classmethod
    def from_xml_path(cls, path: str | PathLike[str]):
        with open(path) as infile:
            return serializer.deserialize_simulation(infile)

    @classmethod
    def from_simulation(cls, simulation: Simulation):
        force = add_imd_force_to_system(simulation.system)
        simulation.context.reinitialize(preserveState=True)
        checkpoint = simulation.context.createCheckpoint()

        pbc_vectors = None
        uses_pbc_wrapping = simulation.system.usesPeriodicBoundaryConditions()
        if uses_pbc_wrapping:
            pbc_vectors = np.array(
                [
                    vector.value_in_unit(nanometer)
                    for vector in simulation.system.getDefaultPeriodicBoxVectors()
                ]
            )

        return LoadedSetup(
            simulation=simulation,
            imd_force=force,
            checkpoint=checkpoint,
            pbc_vectors=pbc_vectors,
            dof=compute_dof(simulation.system),
        )


class Loader[T]:
    setup: LoadedSetup | None = None
    _load: Callable[[T], LoadedSetup]

    def __init__(self, load: Callable[[T], LoadedSetup]):
        self._load = load

    def load(self, config: T) -> LoadedSetup:
        if self.setup is not None:
            return self.setup

        setup = self._load(config)
        self.setup = setup
        return setup
