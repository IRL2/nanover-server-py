import math

import MDAnalysis as mda

from nanover.core import AppServer, Simulation
from nanover.trajectory import FrameData
from nanover.mdanalysis.converter import mdanalysis_to_frame_data


class UniverseSimulation(Simulation):
    @classmethod
    def from_universe(
        cls,
        universe: mda.Universe,
        *,
        name: str | None = None,
    ):
        name = name or universe.filename

        return cls(name=name, universe=universe)

    def __init__(
        self,
        *,
        name: str,
        universe: mda.Universe,
    ):
        self.name = name
        self.universe = universe
        self.app_server: AppServer | None = None

        self.simulation_time = 0
        """Universe time elapsed in playback"""
        self.playback_factor = 1
        """How much Universe time elapses per realtime second"""

    def load(self) -> None:
        pass

    def reset(self, app_server: AppServer) -> None:
        self.app_server = app_server

        self.seek_to_time(0)
        self.app_server.frame_publisher.send_clear()
        self.app_server.frame_publisher.send_frame(self.make_topology_frame())

    def advance_by_one_step(self) -> None:
        self.seek_to_next_frame()

    def advance_by_seconds(self, dt: float) -> None:
        self.seek_to_time(self.simulation_time + dt * self.playback_factor)

    def seek_to_next_frame(self) -> None:
        """Advance simulation time to the time of the next frame and publish it."""
        assert self.app_server is not None

        frame_count = len(self.universe.trajectory)
        frame_length = self.universe.trajectory.dt
        duration = frame_count * frame_length

        # determine previous frame from previous time then take the subsequent frame as next frame and time
        prev_frame = math.floor(frame_count * self.simulation_time / duration)
        next_frame = (prev_frame + 1) % frame_count
        self.simulation_time = next_frame * frame_length

        # update universe frame and publish
        _ = self.universe.trajectory[next_frame]
        self.app_server.frame_publisher.send_frame(self.make_regular_frame())

    def seek_to_time(self, time: float) -> None:
        """Advance simulation time to a specific time and publish the corresponding frame."""
        assert self.app_server is not None

        frame_count = len(self.universe.trajectory)
        frame_length = self.universe.trajectory.dt
        duration = frame_count * frame_length

        # update simulation time and determine corresponding frame
        self.simulation_time = time % duration
        next_frame = math.floor(self.simulation_time / duration * frame_count)

        # update universe frame and publish
        _ = self.universe.trajectory[next_frame]
        self.app_server.frame_publisher.send_frame(self.make_regular_frame())

    def make_topology_frame(self) -> FrameData:
        return mdanalysis_to_frame_data(self.universe, topology=True, positions=True)

    def make_regular_frame(self) -> FrameData:
        return mdanalysis_to_frame_data(self.universe, topology=False, positions=True)
