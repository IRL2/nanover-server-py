import numpy as np
import numpy.typing as npt

from nanover.imd import ParticleInteraction
from nanover.trajectory import FrameData
from nanover.jupyter import ImdAgent
from keyframes import KeyFrame


def fit_keyframe_to_frame(keyframe: KeyFrame, frame: FrameData):
    current = np.array(
        [
            np.mean(frame.particle_positions[target.particles], axis=0)
            for target in keyframe.targets
        ]
    )

    targets = fit_template_points_to_observed(keyframe.centroids, current)

    return current, targets


class SmearAgent(ImdAgent):
    speed = 0.1
    keyframe: KeyFrame | None = None

    def update_interactions(self, full_frame: FrameData, frame_update: FrameData):
        if self.keyframe is None:
            return

        # fit keyframe targets to actual positions
        prev_centroids, next_centroids = fit_keyframe_to_frame(
            self.keyframe, full_frame
        )

        # find necessary motions
        deltas = next_centroids - prev_centroids

        # cap motions by speed
        lengths = np.linalg.norm(deltas, axis=1).reshape(-1, 1)
        cappeds = deltas / lengths
        np.clip(lengths, max=self.speed, out=lengths)
        cappeds *= lengths

        # determine final interaction positions
        target_centroids = prev_centroids + cappeds

        for i, target in enumerate(self.keyframe.targets):
            if lengths[i] < 0.0001:
                continue

            interaction = ParticleInteraction(
                particles=target.particles,
                position=list(target_centroids[i]),
                interaction_type="spring",
                scale=1000,
                max_force=50,
            )
            self.update_interaction(f"interaction.REPLAYER.{i}", interaction)


def fit_template_points_to_observed(
    template_points: npt.NDArray,
    observed_points: npt.NDArray,
):
    # https://web.stanford.edu/class/cs273/refs/umeyama.pdf
    # https://doi.org/10.1109/34.88573
    p = template_points.T
    q = observed_points.T

    # centroids
    cen_P = np.mean(p, axis=1).reshape(-1, 1)
    cen_Q = np.mean(q, axis=1).reshape(-1, 1)

    # centered vectors
    X = p - cen_P
    Y = q - cen_Q

    # svd
    U, sigma, Vt = np.linalg.svd(X @ Y.T)

    # reflection correction
    d = np.identity(Vt.T.shape[1])
    d[-1, -1] = np.linalg.det(Vt.T @ U.T)

    # rotation and translation
    R = Vt.T @ d @ U.T
    t = cen_Q - R @ cen_P

    return (R @ p + t).T
