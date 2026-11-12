from nanover.imd import ParticleInteraction
from nanover.jupyter import ImdAgent
from nanover.trajectory import FrameData
from nanover.utilities.intersection import point_at_distance_along_polyline


class PathFollowerAgent(ImdAgent):
    speed = 0.001
    force_scale = 1000.0
    force_max = 10000.0

    particles: list[int] = []
    path = []
    distance = 0

    def set_particles(self, particles: list[int]):
        self.particles = particles

    def set_path(self, path: list):
        self.path = path
        self.distance = 0

    def update_interactions(self, full_frame: FrameData, frame_update: FrameData):
        result = point_at_distance_along_polyline(self.path, self.distance)

        if result is not None:
            self.interactions.update_interaction(
                "follower",
                ParticleInteraction(
                    position=result.point,
                    particles=self.particles,
                    interaction_type="spring",
                    scale=self.force_scale,
                    force_max=self.force_max,
                ),
            )
            self.distance += self.speed
        else:
            self.interactions.clear()
