from os import PathLike
from pathlib import Path

from typing import Iterable

import MDAnalysis as mda

from nanover.core import AppServer, Simulation
from nanover.trajectory import FrameData
from nanover.mdanalysis.converter import mdanalysis_to_frame_data


class UniverseSimulation(Simulation):
    @classmethod
    def from_path(
        cls,
        structure: PathLike[str],
        *,
        coordinates: Iterable[PathLike[str]] | PathLike[str] | None = None,
        name: str | None = None,
    ):
        """
        Creates UniverseSimulation from path(s).

        :param structure: Path to structure file, such as .pdb or topology.
        :param coordinates: Optional path to trajectory coordinates, if providing a topology to `structure`, this is required.
        :param name: Name of simulation.
        """
        structure = Path(structure)
        name = name or structure.stem

        if coordinates is not None:
            if not isinstance(coordinates, Iterable):
                coordinates = [coordinates]
            elif isinstance(coordinates, str):
                coordinates = Path(coordinates)

            universe = mda.Universe(structure, coordinates)
        else:
            universe = mda.Universe(structure)

        return cls(name=name, universe=universe)

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

    def load(self):
        # Reset simulation to the beginning.
        self._universe_iterator = iter(self.universe)

    def reset(self, app_server: AppServer):
        self.app_server = app_server
        self._next_frame()

        frame = self.make_topology_frame()

        self.app_server.frame_publisher.send_clear()
        self.app_server.frame_publisher.send_frame(frame)

    def _next_frame(self, reset: bool = True) -> None:
        """Advances the internal `Universe` to the next frame or resets to the first timestep if `reset`."""
        try:
            next(self._universe_iterator)
        except StopIteration:
            if reset:
                self._universe_iterator = iter(self.universe)

    def make_topology_frame(self) -> FrameData:
        return mdanalysis_to_frame_data(self.universe, topology=True, positions=True)

    def make_regular_frame(self) -> FrameData:
        return mdanalysis_to_frame_data(self.universe, topology=False, positions=True)

    def advance_by_one_step(self):
        return self.advance_to_next_report()

    def advance_by_seconds(self, dt: float):
        return self.advance_to_next_report()

    def advance_to_next_report(self) -> None:
        assert self.app_server is not None

        self._next_frame()

        # Create frame.
        frame_data = self.make_regular_frame()

        # Send frame
        self.app_server.frame_publisher.send_frame(frame_data)
