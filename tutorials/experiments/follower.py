import numpy as np

from nanover.imd import ParticleInteraction
from nanover.jupyter import ImdAgent
from nanover.trajectory import FrameData


class Path:
    @classmethod
    def empty(cls):
        return cls([])

    @classmethod
    def from_positions(cls, positions: list):
        return cls(positions)

    def __init__(self, positions: list):
        self.positions = positions
        self.lengths = (
            np.linalg.norm(np.diff(positions, axis=0), axis=1)
            if len(positions) >= 2
            else []
        )

    def append_position(self, position: list):
        if self.positions:
            length = np.linalg.norm(np.subtract(self.positions[-1], position))
            self.lengths.append(length)
        self.positions.append(position)

    def position_at_distance(self, distance: float):
        if len(self.positions) == 0:
            return None
        elif len(self.positions) == 1 and distance == 0:
            return self.positions[0]

        target = distance
        for i, length in enumerate(self.lengths):
            if length > target:
                a, b = self.positions[i], self.positions[i + 1]
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
    path = Path.empty()
    distance = 0

    def set_particles(self, particles: list[int]):
        self.particles = particles

    def set_path(self, path: Path):
        self.path = path
        self.distance = 0

    def update_interactions(self, full_frame: FrameData, frame_update: FrameData):
        target = self.path.position_at_distance(self.distance)

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
