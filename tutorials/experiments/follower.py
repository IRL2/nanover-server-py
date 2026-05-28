import numpy as np
import numpy.typing as npt

from nanover.app import OmniRunner
from nanover.imd import ParticleInteraction


from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor
from nanover.utilities.cli import CancellationToken


@dataclass(kw_only=True, eq=False)
class Pin:
    particles: list[int] = field(default_factory=list)
    position: npt.NDArray[np.float32]


@dataclass(kw_only=True, eq=False)
class Checkpoint:
    pins: list[Pin] = field(default_factory=list)


class Path:
    @classmethod
    def from_points(cls, points: list[npt.NDArray[np.float32]]):
        path = cls()
        path.points = points
        path.lengths = np.linalg.norm(np.diff(np.array(points), axis=0), axis=1)
        return path

    def __init__(self):
        self.points: list[npt.NDArray[np.float32]] = []
        self.lengths: list[float] = []

    def add_point(self, point: npt.NDArray[np.float32]):
        if self.points:
            self.lengths.append(np.linalg.norm(point - self.points[-1], axis=0))
        self.points.append(point)

    def get_point_at_distance(self, distance: float) -> npt.NDArray[np.float32] | None:
        for i, length in enumerate(self.lengths):
            if length > distance:
                u = distance / length
                p = self.points[i] + (self.points[i + 1] - self.points[i]) * u
                return p
            distance -= length
        return None


@dataclass(kw_only=True)
class Group:
    particles: list[int] = field(default_factory=list)
    path: Path = field(default_factory=Path)


class Follower:
    @classmethod
    def from_runner(cls, runner: OmniRunner):
        return cls(runner)

    def __init__(self, runner: OmniRunner):
        self._runner = runner
        self._threads = ThreadPoolExecutor(max_workers=1)
        self._cancellation = CancellationToken()
        self._task = None
        self._interactions: set[str] = set()

    def start(self, groups: list[Group], speed: float, release=True):
        if self._task is not None:
            return

        publisher = self._runner.app_server.frame_publisher
        imd = self._runner.app_server.imd
        stream = publisher.subscribe_latest_frames(
            frame_interval=0, cancellation=self._cancellation
        )

        def run():
            distance = 0

            for frame in stream:
                distance += speed

                for i, group in enumerate(groups):
                    key = f"interaction.REPLAYER.{i}"
                    self._interactions.add(key)
                    # centroid = np.average(frame.particle_positions[group.particles], axis=0)

                    target = group.path.get_point_at_distance(distance)

                    if target is None and not release:
                        target = group.path.points[-1]

                    if target is not None:
                        imd.insert_interaction(
                            key,
                            ParticleInteraction(
                                particles=group.particles,
                                position=list(target),
                                interaction_type="spring",
                                scale=500,
                            ),
                        )
                    else:
                        imd.remove_interaction(key)

        self._task = self._threads.submit(run)

    def stop(self):
        self._cancellation.cancel()
        self._threads.shutdown(wait=True)
        for key in self._interactions:
            self._runner.app_server.imd.remove_interaction(key)
