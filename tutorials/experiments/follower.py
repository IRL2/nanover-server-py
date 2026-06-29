import numpy as np

from nanover.imd import ParticleInteraction
from nanover.jupyter import ImdAgent
from nanover.trajectory import FrameData


class Path:
    @classmethod
    def from_points(cls, points: list):
        return cls(points)

    def __init__(self, points: list):
        assert len(points) >= 2
        self._points = points
        self._lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
        self._length = np.sum(self._lengths)

    def point_at_distance(self, distance: float):
        target = distance
        for i, length in enumerate(self._lengths):
            if length > target:
                a, b = self._points[i], self._points[i + 1]
                delta = np.subtract(b, a)
                direction = delta / np.linalg.norm(delta)
                return np.add(a, direction * target)
            target -= length
        else:
            return None


class PathFollowerAgent(ImdAgent):
    speed = 0.01
    force_scale = 100.0
    force_max = 10000.0

    particles: list[int] = []
    path: Path = Path.from_points([[0, 0, 0], [0, 0, 0]])
    distance = 0

    def set_particles(self, particles: list[int]):
        self.particles = particles

    def set_path(self, path: list):
        self.path = Path.from_points(path)
        self.distance = 0

    def update_interactions(self, full_frame: FrameData, frame_update: FrameData):
        target = self.path.point_at_distance(self.distance)

        if target is not None:
            self.interactions.update_interaction(
                "follower",
                ParticleInteraction(
                    position=target,
                    particles=self.particles,
                    interaction_type="spring",
                    scale=self.force_scale,
                    force_max=self.force_max,
                ),
            )
            self.distance += self.speed
        else:
            self.interactions.clear()
