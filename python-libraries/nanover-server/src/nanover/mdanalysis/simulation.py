from os import PathLike
from pathlib import Path

import MDAnalysis as mda

from nanover.core import AppServer, Simulation
from .converter import mdanalysis_to_frame_data


class UniverseSimulation(Simulation):
    @classmethod
    def from_path(
        cls,
        path: PathLike[str],
        *,
        name: str | None = None,
    ):
        path = Path(path)
        name = name or path.stem
        universe = mda.Universe(path)

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

    def load(self): ...

    def seek(self, index: int):
        _ = self.universe.trajectory[index]
        frame = mdanalysis_to_frame_data(self.universe, topology=False)
        self.app_server.frame_publisher.send_frame(frame)

    def reset(self, app_server: AppServer):
        self.app_server = app_server

        frame = mdanalysis_to_frame_data(self.universe)

        self.app_server.frame_publisher.send_clear()
        self.app_server.frame_publisher.send_frame(frame)

    def advance_by_one_step(self): ...

    def advance_by_seconds(self, dt: float): ...
