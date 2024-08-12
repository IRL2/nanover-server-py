from os import PathLike
from pathlib import Path
from typing import Optional, Any

from openmm.app import Simulation, StateDataReporter
from openmm.unit import kilojoule_per_mole

from nanover.app import NanoverImdApplication
from nanover.openmm import serializer, openmm_to_frame_data
from nanover.openmm.imd import (
    create_imd_force,
    NanoverImdReporter,
    add_imd_force_to_system, ImdForceManager, OTHER_FORCE_GROUP_MASK,
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
        self.reporter: Optional[NanoverImdReporter] = None
        self.verbose_reporter: Optional[StateDataReporter] = None

        self.frame_index = 0
        self.imd_force_manager: Optional[ImdForceManager] = None

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
        assert self.simulation is not None and self.checkpoint is not None

        self.app_server = app_server
        self.simulation.context.loadCheckpoint(self.checkpoint)

        try:
            if self.verbose_reporter is not None:
                self.simulation.reporters.remove(self.verbose_reporter)
        except ValueError:
            pass

        if self.verbose_reporter is not None:
            self.simulation.reporters.append(self.verbose_reporter)

        self.imd_force_manager = ImdForceManager(self.app_server.imd, self.imd_force)

        self.app_server.frame_publisher.send_frame(0, self.make_topology_frame(self.simulation))
        self.frame_index = 1

    def advance_by_seconds(self, dt: float):
        self.advance_to_next_report()

    def advance_by_one_step(self):
        self.advance_to_next_report()

    def advance_to_next_report(self):
        assert self.simulation is not None
        self.simulation.step(self.frame_interval)

        step = self.simulation.currentStep
        state = self.simulation.context.getState(
            getPositions=True,
            getForces=self.include_forces,
            getVelocities=self.include_velocities,
            getEnergy=True,
        )

        do_frame = step % self.frame_interval == 0
        do_imd = step % self.force_interval == 0

        positions = state.getPositions(asNumpy=True)
        if do_frame:
            frame_data = self.make_regular_frame(self.simulation, state, positions)
            self.app_server.frame_publisher.send_frame(self.frame_index, frame_data)
            self.frame_index += 1
        if do_imd:
            self.imd_force_manager.update_interactions(self.simulation, positions)

    def make_topology_frame(self, simulation: Simulation):
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        topology = simulation.topology
        frame_data = openmm_to_frame_data(state=state, topology=topology)
        return frame_data

    def make_regular_frame(
        self,
        simulation: Simulation,
        state: State,
        positions: Array2Dfloat,
    ):
        frame_data = openmm_to_frame_data(
            state=state,
            topology=None,
            include_positions=False,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
        )
        frame_data.particle_positions = positions
        self.imd_force_manager.add_to_frame_data(frame_data)

        # Get the simulation state excluding the IMD force, and recalculate potential energy without it:
        energy_no_imd = (
            simulation.context.getState(getEnergy=True, groups=OTHER_FORCE_GROUP_MASK)
            .getPotentialEnergy()
            .value_in_unit(kilojoule_per_mole)
        )
        frame_data.potential_energy = energy_no_imd

        return frame_data
