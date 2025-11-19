from abc import abstractmethod
from typing import Protocol, TypeVar, Callable

from . import AppServer
from nanover.trajectory import FrameData


S = TypeVar("S")


class Simulation(Protocol[S]):
    name: str

    def load(self): ...
    def reset(self, app_server: AppServer): ...  # Change to `reset(self)` ?

    @abstractmethod
    def make_topology_frame(self) -> FrameData: ...

    @abstractmethod
    def make_regular_frame(self) -> FrameData: ...

    @abstractmethod
    def advance_by_one_step(self):
        # Should modify behaviour of `_advance_to_next_report` instead i.e. so moves by steps than sec.
        ...

    @abstractmethod
    def advance_by_seconds(self, dt: float):
        # Should modify behaviour of `_advance_to_next_report` instead i.e. so moves by sec rather than step.
        ...

    @abstractmethod
    def _advance_to_next_report(self) -> S:
        # Unified method to get next relevant simulation object
        # e.g. return self.simulation after .step() for openmm simulation
        ...

    @abstractmethod
    def _frame_from_simulation(self, current_timepoint: S) -> FrameData:
        """Constructs a frame from the current simulation timepoint."""
        ...

    @abstractmethod
    def _send_frame(self, frame: FrameData) -> None:
        # Sends frame via frame_publisher. Requred to maintain no AppServer member.
        ...

    @abstractmethod
    def _cleanup(self, most_recent: FrameData | None = None) -> None:
        """Cleans up any parameters/settings used in stepping to the next timepoint."""
        ...

    def advance_to_next_report(self) -> None:
        """Advance the underlying simulation to the next frame and send it over the network."""
        current_timestep = self._advance_to_next_report()
        current_frame = self._frame_from_simulation(current_timestep)

        self._send_frame(current_frame)
        self._cleanup(current_frame)


class ReportableSimulation(Simulation):
    _measurements = list[float]
    measure_callbacks: set[Callable[[FrameData], float] | Callable[[S], float]]

    def add_measurement(
        self, callback: Callable[[S], float] | Callable[[FrameData], float]
    ) -> None: ...

    def _process_measurements(self, current_state: S) -> None:
        _measurements = [func(current_state) for func in self.measure_callbacks]
        # Add new measurement data to frame
        # _measurements.add_to_frame(frame)

    def advance_to_next_report(self):
        current_timestep = self._advance_to_next_report()
        current_frame = self._frame_from_simulation(current_timestep)
        # self._process_measurements(current_timestep).add_to_frame(current_frame)

        self._send_frame(current_frame)
        self._cleanup(current_frame)
