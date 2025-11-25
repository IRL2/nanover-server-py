import MDAnalysis as mda

from nanover.core import AppServer, Simulation
from nanover.trajectory import FrameData
from nanover.mdanalysis.converter import mdanalysis_to_frame_data


class StepTracker:
    """Class to store step sizes for Simulation class advancement."""

    def __init__(self, dt: float | None = None):
        self.default_dt: float | None = dt
        "Target time to step in each advancement, `None` implies stepping 1 frame."
        self.current_dt: float | None = None
        "Target time to step for the next advancement only, `None` implies stepping 1 frame."

    def next_ts(self) -> float | None:
        """Returns the next target time step or `None` if stepping 1 frame."""
        return self.current_dt or self.default_dt

    def __next__(self) -> float | None:
        """Yields next target time step."""
        next_ts = self.default_dt

        # Exhaust current_dt if it exists.
        if self.current_dt is not None:
            next_ts = self.current_dt
            self.current_dt = None
        return next_ts


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
        self._universe_iterator = iter(self.universe.trajectory)
        self._time_to_step: StepTracker = StepTracker()
        self.time_mapping_factor = (
            1 / 3
        )  # factor assumes rw and sim time are s and ps, respectively.
        "Factor to map time elapsed in real world (s) and simulation (ps) when advancing by a given time."

        self.app_server: AppServer | None = None

    def load(self):
        # Reset simulation iterator to the beginning.
        self._universe_iterator = iter(self.universe.trajectory)

    def reset(self, app_server: AppServer):
        """Reset the simulation to the first timestep then sends a frame."""
        self.app_server = app_server
        self._universe_iterator = iter(self.universe.trajectory)
        next(self._universe_iterator)  # Advance to first frame.

        frame = self.make_topology_frame()

        self.app_server.frame_publisher.send_clear()
        self.app_server.frame_publisher.send_frame(frame)

    def _next_frame(self, reset: bool = True) -> None:
        """
        Advances the internal `Universe` to the next frame.
        Resets to the first timestep if finished iterating over the simulation and `reset`.
        """
        try:
            # Simply advance to next step if not advancing by a given amount of time.
            if (next_ts := next(self._time_to_step)) is None:
                next(self._universe_iterator)
            else:
                # Get number of frames to advance - may slightly undershoot if not exact multiple of `dt`.
                num_frames_to_advance = int(
                    (next_ts / self.time_mapping_factor) // self._universe_iterator.dt
                )
                for _ in range(num_frames_to_advance):
                    next(self._universe_iterator)
        except StopIteration:
            if reset:
                self._universe_iterator = iter(self.universe.trajectory)
                next(self._universe_iterator)  # Advance to first frame.

    def make_topology_frame(self) -> FrameData:
        return mdanalysis_to_frame_data(self.universe, topology=True, positions=True)

    def make_regular_frame(self) -> FrameData:
        return mdanalysis_to_frame_data(self.universe, topology=False, positions=True)

    @property
    def time_to_step(self) -> float | None:
        """
        The default amount of time (ps) the simulation will attempt to step when advancing.
        Returns `None` if there is no target time, indicating it will advance to the next timestep instead.
        """
        return self._time_to_step.default_dt

    @time_to_step.setter
    def time_to_step(self, dt: float) -> None:
        self._time_to_step.default_dt = dt

    def advance_by_one_step(self) -> None:
        """Advance the simulation to the next timestep."""
        self._time_to_step.current_dt = None
        return self.advance_to_next_report()

    def advance_by_seconds(self, dt: float) -> None:
        """
        Advance the simulation by a given time change `dt` in real world seconds.
        This is mapped to a simulation time change with `time_mapping_factor` member.
        """
        self._time_to_step.current_dt = dt
        return self.advance_to_next_report()

    def advance_to_next_report(self) -> None:
        """
        Advances the simulation to the next desired time point and sends the new state as a frame.
        If `advance_by_seconds` has been called, will attempt to progress to the timestep
        corresponding to the given change in time. Otherwise sends the next timepoint.
        """
        assert self.app_server is not None

        self._next_frame()

        # Create frame.
        frame_data = self.make_regular_frame()

        # Send frame
        self.app_server.frame_publisher.send_frame(frame_data)
