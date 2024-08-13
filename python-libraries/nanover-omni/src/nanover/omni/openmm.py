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
    NON_IMD_FORCES_GROUP_MASK,
)
from nanover.trajectory.frame_data import Array2Dfloat


class OpenMMSimulation:
    """
    A wrapper for OpenMM simulations so they can be run inside the OmniRunner.
    """

    @classmethod
    def from_simulation(cls, simulation: Simulation, *, name: Optional[str] = None):
        """
        Construct this from an existing OpenMM simulation.
        :param simulation: An existing OpenMM Simulation
        :param name: An optional name for the simulation
        :return:
        """
        sim = cls(name)
        sim.simulation = simulation
        sim.imd_force = add_imd_force_to_system(simulation.system)
        sim.simulation.context.reinitialize(preserveState=True)

        sim.checkpoint = sim.simulation.context.createCheckpoint()

        return sim

    @classmethod
    def from_xml_path(cls, path: PathLike[str], *, name: Optional[str] = None):
        """
        Construct this from an existing NanoVer OpenMM XML file at a given path.
        :param path: Path of the NanoVer OpenMM XML file
        :param name: An optional name for the simulation
        :return:
        """
        sim = cls(name or Path(path).stem)
        sim.xml_path = path
        return sim

    def __init__(self, name: Optional[str] = None):
        self.name = name or "Unnamed OpenMM Simulation"

        self.xml_path: Optional[PathLike[str]] = None
        self.app_server: Optional[NanoverImdApplication] = None

        self.frame_interval = 5
        self.include_velocities = False
        self.include_forces = False
        self.platform: Optional[str] = None

        self.imd_force = create_imd_force()
        self.simulation: Optional[Simulation] = None
        self.checkpoint: Optional[Any] = None
        self.verbose_reporter: Optional[StateDataReporter] = None

        self.frame_index = 0
        self.imd_force_manager: Optional[ImdForceManager] = None

    def load(self):
        """
        Load and set up the simulation if it isn't done already.
        """
        if self.xml_path is None or self.simulation is not None:
            return

        with open(self.xml_path) as infile:
            self.imd_force = create_imd_force()
            self.simulation = serializer.deserialize_simulation(
                infile, imd_force=self.imd_force, platform_name=self.platform
            )

        self.checkpoint = self.simulation.context.createCheckpoint()

    def reset(self, app_server: NanoverImdApplication):
        """
        Reset the simulation to its initial conditions, reset IMD interactions, and reset frames to begin with topology
        and continue.
        :param app_server: The app server hosting the frame publisher and imd state
        :return:
        """
        assert self.simulation is not None and self.checkpoint is not None

        self.app_server = app_server
        self.simulation.context.loadCheckpoint(self.checkpoint)
        self.imd_force_manager = ImdForceManager(self.app_server.imd, self.imd_force)

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

    def advance_by_seconds(self, dt: float):
        """
        Advance the simulation playback by some seconds. This time is ignored and the next frame is generated and sent.
        :param dt: Time to advance playback by
        """
        self.advance_to_next_report()

    def advance_by_one_step(self):
        """
        Advance the simulation to the next point a frame should be reported, and send that frame.
        :return:
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
            enforcePeriodicBox=False,
        )
        positions = state.getPositions(asNumpy=True)

        # update imd forces and energies
        self.imd_force_manager.update_interactions(self.simulation, positions)

        # generate the next frame with the existing (still valid) positions
        frame_data = self.make_regular_frame(positions)

        # send the next frame
        self.app_server.frame_publisher.send_frame(self.frame_index, frame_data)
        self.frame_index += 1

    def make_topology_frame(self):
        """
        Make a NanoVer FrameData corresponding to the current particle positions and topology of the simulation.
        """
        assert self.simulation is not None

        state = self.simulation.context.getState(getPositions=True, getEnergy=True)
        topology = self.simulation.topology
        frame_data = openmm_to_frame_data(state=state, topology=topology)
        return frame_data

    def make_regular_frame(self, positions: Optional[Array2Dfloat] = None):
        """
        Make a NanoVer FrameData corresponding to the current state of the simulation.
        :param positions: Optionally provided particle positions to save fetching them again.
        """
        assert self.simulation is not None and self.imd_force_manager is not None

        # fetch omm state
        state = self.simulation.context.getState(
            getPositions=positions is None,
            getForces=self.include_forces,
            getVelocities=self.include_velocities,
            getEnergy=True,
            groups=NON_IMD_FORCES_GROUP_MASK,
        )

        # generate frame based on basic omm state
        frame_data = openmm_to_frame_data(
            state=state,
            topology=None,
            include_positions=positions is None,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
        )

        # add any provided positions
        if positions is not None:
            frame_data.particle_positions = positions

        # add imd force information
        self.imd_force_manager.add_to_frame_data(frame_data)

        return frame_data
