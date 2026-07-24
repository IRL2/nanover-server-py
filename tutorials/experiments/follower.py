import numpy as np
import numpy.typing as npt

from dataclasses import dataclass, field

from nanover.imd import ParticleInteraction
from nanover.trajectory import FrameData
from nanover.jupyter import ImdAgent


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


class Follower(ImdAgent):
    groups: list[Group] | None = None
    speed = 0.1
    release = True
    distance = 0

    def update_interactions(self, full_frame: FrameData, frame_update: FrameData):
        if self.groups is None:
            return

        self.distance += self.speed

        for i, group in enumerate(self.groups):
            key = f"interaction.REPLAYER.{i}"
            target = group.path.get_point_at_distance(self.distance)

            if target is None and not self.release:
                target = group.path.points[-1]

            if target is not None:
                self.update_interaction(
                    key,
                    ParticleInteraction(
                        particles=group.particles,
                        position=list(target),
                        interaction_type="spring",
                        scale=500,
                    ),
                )
            else:
                self.remove_interaction(key)
