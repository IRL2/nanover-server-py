import math

import MDAnalysis as mda
from nanover.core import AppServer, Simulation
from nanover.mdanalysis.converter import mdanalysis_to_frame_data
from nanover.trajectory import FrameData


class UniverseSimulation(Simulation):
    @classmethod
    def from_universe(
        cls,
        universe: mda.Universe,
        *,
        name: str | None = None,
    ):
        name = name or str(universe.filename)

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

        self.simulation_time = 0
        _ = self.universe.trajectory[0]

        self.app_server.frame_publisher.send_clear()
        self.app_server.frame_publisher.send_frame(self.make_topology_frame())

    def advance_by_one_step(self) -> None:
        if self.frame_count > 1:
            self._seek_to_next_frame()

    def advance_by_seconds(self, dt: float) -> None:
        if self.frame_count > 1:
            self._seek_to_time(self.simulation_time + dt * self.playback_factor)

    @property
    def frame_count(self):
        return len(self.universe.trajectory)

    @property
    def frame_length(self):
        return self.universe.trajectory.dt

    @property
    def duration(self):
        return self.frame_count * self.frame_length

    @property
    def prev_frame(self):
        return math.floor(self.frame_count * self.simulation_time / self.duration)

    def seek_to_frame_index(self, index: int) -> None:
        assert self.app_server is not None

        next_frame = index

        # align time to start of target frame
        self.simulation_time = next_frame * self.frame_length

        # update universe frame and publish
        _ = self.universe.trajectory[next_frame]
        self.app_server.frame_publisher.send_frame(self.make_regular_frame())

    def _seek_to_next_frame(self) -> None:
        """Advance simulation time to the time of the next frame and publish it."""
        assert self.app_server is not None

        # determine next frame index then seek to it
        next_frame = (self.prev_frame + 1) % self.frame_count
        self.seek_to_frame_index(next_frame)

    def _seek_to_time(self, time: float) -> None:
        """Advance simulation time to a specific time and publish the corresponding frame."""
        assert self.app_server is not None

        # update simulation time and determine corresponding frame
        self.simulation_time = time % self.duration
        next_frame = math.floor(self.simulation_time / self.duration * self.frame_count)

        # update universe frame and publish
        _ = self.universe.trajectory[next_frame]
        self.app_server.frame_publisher.send_frame(self.make_regular_frame())

    def make_topology_frame(self) -> FrameData:
        return mdanalysis_to_frame_data(self.universe, topology=True, positions=True)

    def make_regular_frame(self) -> FrameData:
        return mdanalysis_to_frame_data(self.universe, topology=False, positions=True)
