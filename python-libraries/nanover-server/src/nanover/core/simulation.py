from abc import abstractmethod
from typing import Protocol, TypeVar

from . import AppServer
from nanover.trajectory import FrameData


S = TypeVar("S")


class Simulation(Protocol[S]):
    name: str
    app_server: AppServer

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
        # e.g. return self.simulation.step() for openmm simulation
        ...

    @abstractmethod
    def _frame_from_simulation(self, current_timepoint: S) -> FrameData:
        """Constructs a frame from the current simulation timepoint."""
        ...

    @abstractmethod
    def _cleanup(self, most_recent: FrameData | None = None) -> None:
        """Cleans up any parameters/settings used in stepping to the next timepoint."""
        ...

    def send_next_frame(self) -> None:
        """Advance the underlying simulation to the next frame and send it over the network."""
        current_timestep = self._advance_to_next_report()
        current_frame = self._frame_from_simulation(current_timestep)

        self.app_server.frame_publisher.send_frame(current_frame)
        self._cleanup(current_frame)
