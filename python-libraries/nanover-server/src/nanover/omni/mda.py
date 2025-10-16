from os import PathLike
from pathlib import Path

import MDAnalysis as mda

from nanover.app.types import AppServer
from nanover.mdanalysis import mdanalysis_to_frame_data


class UniverseSimulation:
    @classmethod
    def from_path(
        cls,
        *,
        name: str | None = None,
        path: PathLike[str],
    ):
        path = Path(path)
        name = name or path.stem
        universe = mda.Universe(path)

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

    def reset(self, app_server: AppServer):
        self.app_server = app_server

        frame = mdanalysis_to_frame_data(self.universe)

        self.app_server.frame_publisher.send_clear()
        self.app_server.frame_publisher.send_frame(frame)

    def advance_by_one_step(self): ...

    def advance_by_seconds(self, dt: float): ...
