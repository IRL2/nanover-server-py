from os import PathLike
from pathlib import Path
from typing import Optional, Any

from openmm.app import Simulation, StateDataReporter

from nanover.app import NanoverImdApplication
from nanover.openmm import serializer, openmm_to_frame_data
from nanover.openmm.imd import (
    create_imd_force,
    add_imd_force_to_system,
    ImdForceManager,
)
from nanover.protocol.state import State
from nanover.trajectory.frame_data import Array2Dfloat


class OpenMMSimulation:
    @classmethod
    def from_simulation(cls, simulation: Simulation, *, name: Optional[str] = None):
        sim = cls(name)
        sim.simulation = simulation
        sim.imd_force = add_imd_force_to_system(simulation.system)
        sim.simulation.context.reinitialize(preserveState=True)

        sim.checkpoint = sim.simulation.context.createCheckpoint()

        return sim

    @classmethod
    def from_xml_path(cls, path: PathLike[str], *, name: Optional[str] = None):
        sim = cls(name or Path(path).stem)
        sim.xml_path = path
        return sim

    def __init__(self, name: Optional[str] = None):
        self.name = name or "Unnamed OpenMM Simulation"

        self.xml_path: Optional[PathLike[str]] = None
        self.app_server: Optional[NanoverImdApplication] = None

        self.frame_interval = 5
        self.force_interval = 5
        self.include_velocities = False
        self.include_forces = False
        self.platform: Optional[str] = None

        self.imd_force = create_imd_force()
        self.simulation: Optional[Simulation] = None
        self.checkpoint: Optional[Any] = None
        self.verbose_reporter: Optional[StateDataReporter] = None

        self._frame_index = 0
        self._imd_force_manager: Optional[ImdForceManager] = None

    def load(self):
        if self.xml_path is None or self.simulation is not None:
            return

        with open(self.xml_path) as infile:
            self.imd_force = create_imd_force()
            self.simulation = serializer.deserialize_simulation(
                infile, imd_force=self.imd_force, platform_name=self.platform
            )

        self.checkpoint = self.simulation.context.createCheckpoint()

    def reset(self, app_server: NanoverImdApplication):
        assert (
            self.simulation is not None
            and self.checkpoint is not None
            and self._imd_force_manager is not None
        )

        self.app_server = app_server
        self.simulation.context.loadCheckpoint(self.checkpoint)

        try:
            if self.verbose_reporter is not None:
                self.simulation.reporters.remove(self.verbose_reporter)
        except ValueError:
            pass

        if self.verbose_reporter is not None:
            self.simulation.reporters.append(self.verbose_reporter)

        self._imd_force_manager = ImdForceManager(self.app_server.imd, self.imd_force)
        frame = self.make_topology_frame(self.simulation)
        self.app_server.frame_publisher.send_frame(0, frame)
        self._frame_index = 1

    def make_topology_frame(self, simulation: Simulation):
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        topology = simulation.topology
        frame_data = openmm_to_frame_data(state=state, topology=topology)
        return frame_data

    def make_regular_frame(self, state: State, positions: Array2Dfloat):
        frame_data = openmm_to_frame_data(
            state=state,
            topology=None,
            include_positions=False,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
        )
        frame_data.particle_positions = positions
        self._imd_force_manager.add_to_frame_data(frame_data)

        return frame_data

    def advance_to_next_report(self):
        assert (
            self._imd_force_manager is not None
            and self.simulation is not None
            and self._imd_force_manager is not None
        )
        self.simulation.step(self.frame_interval)

        state = self.simulation.context.getState(
            getPositions=True,
            getVelocities=self.include_velocities,
            getForces=self.include_forces,
            getEnergy=True,
        )
        positions = state.getPositions(asNumpy=True)
        if self.simulation.currentStep % self.frame_interval == 0:
            frame_data = self.make_regular_frame(state, positions)
            self.app_server.frame_publisher.send_frame(self._frame_index, frame_data)
            self._frame_index += 1
        if self.simulation.currentStep % self.force_interval == 0:
            self._imd_force_manager.update_interactions(self.simulation, positions)

    def advance_by_seconds(self, dt: float):
        self.advance_to_next_report()

    def advance_by_one_step(self):
        self.advance_to_next_report()
